#include <kmerinshort.hpp>
#include <stdio.h>

using namespace std;
typedef uint32_t count_t;

/********************************************************************************/

#define Saturated_inc8(a)  ((a == 0xFF ) ? 0xFF : a++)
#define Saturated_inc16(a)  ((a == 0xFFFF ) ? 0xFFFF : a++)
#define Saturated_inc32(a)  ((a == 0xFFFFFFFF ) ? 0xFFFFFFFF : a++)

const char * kis::STR_QUERY = "-query";

#pragma mark -
#pragma mark KmerCounter

template <typename modelT>
class KmerCounter
{
	typedef typename modelT::Iterator KmerIterator;
	
public:
	
	void operator() (Sequence& sequence)
	{
		
		//printf("sequence index %zu \n",sequence.getIndex());
		
		KmerIterator itKmer (*_model);
		itKmer.setData (sequence.getData());
		
		
		for (itKmer.first(); !itKmer.isDone(); itKmer.next())
		{
			
			if(_count_table_current_bank_global [itKmer->value().getVal()]< 4294967295)
				__sync_fetch_and_add (& _count_table_current_bank_global [itKmer->value().getVal()], 1);
			
		}
		
		seqdone++;
		if (seqdone > 1000)   {  _progress.inc (seqdone);  seqdone = 0;  }
		
	}
	
	//constructor
	KmerCounter(modelT * model,
				count_t * count_table_current_bank,
				u_int64_t nbDiffKmers,
				IteratorListener* progress,
				std::string outputFilename)
	:_model(model),seqdone(0),_progress(progress,System::thread().newSynchronizer()),_nbDiffKmers(nbDiffKmers),_count_table_current_bank_global(count_table_current_bank),_outputFilename(outputFilename)
	{
		
	}
	
private:
	
	modelT * _model;
	int seqdone;
	ProgressSynchro  _progress;
	u_int64_t _nbDiffKmers;
	
	count_t * _count_table_current_bank_global;
	std::string _outputFilename;
	
};


#pragma mark -
#pragma mark kmerinshort


kis::kis ()  : Tool ("kis"),_progress (0)
{

	
	getParser()->push_back (new OptionOneParam (STR_URI_FILE, "input file ",   true));
	getParser()->push_back (new OptionOneParam (STR_KMER_SIZE, "ksize",   true));
	getParser()->push_back( new OptionOneParam(STR_URI_OUTPUT, "output file", false));
	getParser()->push_back(new  OptionNoParam("-dont-reverse", "do not reverse kmers, count forward and reverse complement separately", false));
	
	getParser()->push_back( new OptionOneParam(STR_QUERY, "query", false));

	
}

static const size_t span = KMER_SPAN(0);


#define IX(kmer,idx) ((kmer)+(_nbDiffKmers)*(idx))

void kis::query()
{
	//fread ..
	
	fread(&_nbDiffKmers,sizeof(_nbDiffKmers),1,_countfile);
	
	printf("_nbDiffKmers %llu \n",_nbDiffKmers);
	
	count_t * count_table;
 

	count_table = (count_t *) calloc(_nbDiffKmers, sizeof(count_t));
	fread(count_table,_nbBanks * _nbDiffKmers , sizeof(count_t),_countfile);
	
	
	//query input file
	
	_inputBank = Bank::open(_inputFilename);

	//wrapped with a progress iterator
	ProgressIterator<Sequence> *it = new ProgressIterator<Sequence> (*_inputBank, "Iterating sequences");
	
	// We declare a kmer model with a given span size.
	Kmer<span>::ModelCanonical model (_kmerSize);
	
	// We declare a kmer iterator
	Kmer<span>::ModelCanonical::Iterator itKmer (model);
	
	
	// We loop over sequences.
	for (it->first(); !it->isDone(); it->next())
	{
		// We set the data from which we want to extract kmers.
		itKmer.setData ((*it)->getData());
		// We iterate the kmers.
		for (itKmer.first(); !itKmer.isDone(); itKmer.next())
		{
			 cout << model.toString (itKmer->value()) <<  " abundance in first file : "  << count_table[itKmer->value().getVal()] << endl;
		}
	}
	
	
	free(count_table);
}

void kis::execute ()
{
	_query_mode =false;
	_dontreverse = getParser()->saw("-dont-reverse");
	if(_dontreverse)
		fprintf(stderr,"Counting kmers on fasta strand only\n");
	else
		fprintf(stderr,"Counting in canonical mode (kmer and their reverse complement counted together) \n");
	
	
	getInput()->setInt (STR_VERBOSE, 0); //we force verbose 0
	
	
	_inputFilename = getInput()->getStr(STR_URI_FILE);
	
	
	_kmerSize      = getInput()->getInt (STR_KMER_SIZE);
	_inputBank = Bank::open(_inputFilename);
	
	
	if(getParser()->saw(STR_QUERY))
	{
		_query_mode = true;
		_countFilename = getInput()->getStr(STR_QUERY);
		_countfile = fopen(_countFilename.c_str(),"r");
		if(_countfile==NULL)
		{
			fprintf(stderr,"cannot open %s \n",_countFilename.c_str());
			exit(1);
		}
	
		
		query();
		
		return;
	}
	
	
	if(getParser()->saw(STR_URI_OUTPUT))
	{
		_outputFilename = getInput()->getStr(STR_URI_OUTPUT);
		_outfile = fopen(_outputFilename.c_str(),"w");
		
		if(_outfile==NULL)
		{
			fprintf(stderr,"cannot open %s \n",_outputFilename.c_str());
			exit(1);
		}
	}
	else
	{
		_outfile = stdout;
	}
	

	
	u_int64_t number, totalSize, maxSize;
	_inputBank->estimate (number, totalSize, maxSize);
	ProgressTimer * pt =  new ProgressTimer(number, "counting kmers");
	setProgress (new ProgressSynchro (pt,System::thread().newSynchronizer()));
	
	
	ModelCanonical model(_kmerSize);
	ModelDirect modeld(_kmerSize);
	
	
	Iterator<Sequence>* it = _inputBank->iterator();
	std::vector<Iterator<Sequence>*> itBanks =  it->getComposition();
	_nbBanks = itBanks.size();
	_nb_cores = getInput()->getInt(STR_NB_CORES);
	
	
	_nbDiffKmers = model.getKmerMax().getVal() +1 ;
	
	
	count_t * count_table;
 
	
	// count together all input files :
	_nbBanks = 1;
	
	count_table = (count_t *) malloc( _nbBanks * _nbDiffKmers * sizeof(count_t));
	memset(count_table,0, _nbBanks * _nbDiffKmers * sizeof(count_t));
	
	//printf("nb diff kmers  %llu  nb banks %lu  memory allocated %llu MB\n",_nbDiffKmers,itBanks.size(), _nbBanks * _nbDiffKmers * sizeof(count_t) /1024LL/1024LL );
	
	count_t * count_table_current_bank;
	
	
	_progress->init ();
	
	//loop over all input banks
	for (size_t ii=0; ii<itBanks.size(); ii++)
	{
		
		
	//	if(itBanks.size() > 1 )
	//	{
	//		fprintf(_outfile,"-----Bank %zu ------\n",ii);
	//	}
		
		Iterator<Sequence>* itSeq = itBanks[ii];
		
		
		// count together all input files :
		count_table_current_bank = count_table;
		
		//or separately
		//		count_table_current_bank =  & (count_table[ _nbDiffKmers * ii]);
		
		//loop over sequences.
		setDispatcher (  new Dispatcher (_nb_cores) );
		
		if(_dontreverse)
			getDispatcher()->iterate (itSeq,  KmerCounter<ModelDirect>(&modeld,count_table_current_bank,_nbDiffKmers,_progress,_outputFilename), 1000);
		else
			getDispatcher()->iterate (itSeq,  KmerCounter<ModelCanonical>(&model,count_table_current_bank,_nbDiffKmers,_progress,_outputFilename), 1000);
		
		//output result
	}
	
	_progress->finish ();
	
	
	//save in raw binary:
	
	fwrite(&_nbDiffKmers,sizeof(_nbDiffKmers),1,_outfile);
	fwrite(count_table,_nbBanks * _nbDiffKmers , sizeof(count_t),_outfile);
	
	
	//in ascii
	/*
	// count together all input files :
	 for(uint64_t km = 0; km <_nbDiffKmers; km++ )
	 {
		kmer_type km2; km2.setVal(km);
		fprintf(_outfile,"%s",model.toString(km2).c_str());
	 	fprintf(_outfile,"\t%i", count_table[IX(km,0)]);
		fprintf(_outfile,"\n");
	 }
	*/
	
	
	if(getParser()->saw(STR_URI_OUTPUT) )
		fclose(_outfile);
	
	
	free(count_table);
	
}








