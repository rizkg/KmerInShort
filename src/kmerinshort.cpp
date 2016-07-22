#include <kmerinshort.hpp>
using namespace std;
typedef uint32_t count_t;

/********************************************************************************/

#define Saturated_inc8(a)  ((a == 0xFF ) ? 0xFF : a++)
#define Saturated_inc16(a)  ((a == 0xFFFF ) ? 0xFFFF : a++)
#define Saturated_inc32(a)  ((a == 0xFFFFFFFF ) ? 0xFFFFFFFF : a++)

const char* kis::STR_REFK = "-kval";

#pragma mark -
#pragma mark KmerCounter

template <typename modelT>
class KmerCounter
{
typedef typename modelT::Iterator KmerIterator;

public:
	
	void operator() (Sequence& sequence)
	{
		//raz du cpt des kmer
		if(_kismode)
		{
			_count_table_current_bank.assign(_nbDiffKmers,0);
		}
		
		//printf("sequence index %zu \n",sequence.getIndex());
		
		uint64_t nbkmers = 0;
		
		KmerIterator itKmer (*_model);
		itKmer.setData (sequence.getData());

		int pas = _step - 1 ;
		int off = 0;
		
		for (itKmer.first(); !itKmer.isDone(); itKmer.next())
		{
			if(off < _offset)
			{
				off++;
				continue;
			}
			pas ++;
			if(pas != _step)
			{
				continue;
			}
			
			nbkmers ++;
			
			if(_kismode)
			{
				if(_count_table_current_bank [itKmer->value().getVal()]< 4294967295)
				{
					_count_table_current_bank [itKmer->value().getVal()]++;
				}
			}
			else //global cpt, needs to sync
			{
				if(_count_table_current_bank_global [itKmer->value().getVal()]< 4294967295)
					__sync_fetch_and_add (& _count_table_current_bank_global [itKmer->value().getVal()], 1);
				
			}
			
			pas = 0;
		}
		
		if(_kismode)
		{
			
			double sequence_score = 0;
			double kmer_score;
			//now compute the sequence score
			for(uint64_t km = 0; km <_nbDiffKmers; km++ )
			{
				kmer_score = _count_table_current_bank[km]  * _ref_table [km];
				sequence_score += kmer_score;
			}
			
			sequence_score = sequence_score  / nbkmers;
			
			
			//store resu
			_resuvec_perseq [sequence.getIndex() ] =  sequence_score;
			
		}
		
		seqdone++;
		if (seqdone > 1000)   {  _progress.inc (seqdone);  seqdone = 0;  }
		
	}
	
	//constructor
	KmerCounter(modelT * model,
				count_t * count_table_current_bank,
				u_int64_t nbDiffKmers,
				IteratorListener* progress,int offset, int step, double * ref_table, double * resuvec_perseq, bool kismode)
	:_model(model),seqdone(0),_progress(progress,System::thread().newSynchronizer()),_offset(offset),_step(step), _ref_table(ref_table),_resuvec_perseq(resuvec_perseq),_nbDiffKmers(nbDiffKmers),_kismode(kismode),
	_count_table_current_bank_global(count_table_current_bank)
	{
		if(_kismode)
		{
			_count_table_current_bank.resize(_nbDiffKmers,0);
		}
	}
	

	private:

	
	modelT * _model;
	int seqdone;
	ProgressSynchro  _progress;
	int _offset;
	int _step;
	double * _ref_table;
	double * _resuvec_perseq;
	u_int64_t _nbDiffKmers;
	bool _kismode;

	std::vector<count_t> _count_table_current_bank; //now per seq so per thread
	count_t * _count_table_current_bank_global;

	
};


#pragma mark -
#pragma mark kmerinshort


kis::kis ()  : Tool ("kis"),_progress (0)
{
	_offset =0;
	_step = 1;

	getParser()->push_back (new OptionOneParam (STR_URI_FILE, "input file ",   true));
	getParser()->push_back (new OptionOneParam (STR_KMER_SIZE, "ksize",   true));
	getParser()->push_back( new OptionOneParam(STR_URI_OUTPUT, "output file", false));
	getParser()->push_back (new OptionOneParam ("-offset", "starting offset",   false,"0"));
	getParser()->push_back (new OptionOneParam ("-step", "step",   false,"1"));
	getParser()->push_back (new OptionOneParam (STR_REFK, "file with kmer values ",   false));
	getParser()->push_back(new  OptionNoParam("-dont-reverse", "do not reverse kmers, count forward and reverse complement separately", false));

}

static const size_t span = KMER_SPAN(0);


#define IX(kmer,idx) ((kmer)+(_nbDiffKmers)*(idx))


void kis::parse_ref ()
{
	//printf("nb diff kmers %lli \n",_nbDiffKmers);
	
	char *line = NULL;
	size_t linecap = 0;
	ssize_t linelen;
	double new_val;
	u_int64_t cpt = 0;
	
	while ((linelen = getline(&line, &linecap, _kmerRefFile)) > 0)
	{
		if(cpt >= _nbDiffKmers)
		{
			fprintf(stderr,"ERROR: too many lines in ref file, must have %lli (4^%zu) lines \n",_nbDiffKmers,_kmerSize);
			exit(1);
		}
		
		sscanf(line,"%lf",&new_val);
		_ref_table[cpt] = new_val;
		cpt ++;
	}
	
	//printf("cpt %llu / %lli\n",cpt,_nbDiffKmers);
	if(cpt < _nbDiffKmers)
	{
		fprintf(stderr,"ERROR: not enough lines in ref file, must have %lli (4^%zu) lines \n",_nbDiffKmers,_kmerSize);
		exit(1);
	}
}


void kis::execute ()
{
	_dontreverse = getParser()->saw("-dont-reverse");
	if(_dontreverse)
		fprintf(stderr,"Counting kmers on fasta strand only\n");
	else
		fprintf(stderr,"Counting in canonical mode (kmer and their reverse complement counted together) \n");

	 getInput()->setInt (STR_VERBOSE, 0); //we force verbose 0
	
	_kismode = getParser()->saw(STR_REFK);
	
	_offset = (getInput()->getInt("-offset"));
	if(_offset<0)
	{
		fprintf(stderr,"offset should not be <0 \n");
		exit(1);
	}
	_step = (getInput()->getInt("-step"));
	
	if(_step<0)
	{
		fprintf(stderr,"step should not be <0 \n");
		exit(1);
	}
	
	_inputFilename = getInput()->getStr(STR_URI_FILE);
	if(_kismode)
	{
		_kmerValFilename = getInput()->getStr(STR_REFK);
		_kmerRefFile = fopen(_kmerValFilename.c_str(),"r");
		
		if(_kmerRefFile==NULL)
		{
			fprintf(stderr,"cannot open %s \n",_kmerValFilename.c_str());
			exit(1);
		}
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
	
	_kmerSize      = getInput()->getInt (STR_KMER_SIZE);
	_inputBank = Bank::open(_inputFilename);
	
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
	
	//load ref values in ref_table
	if(_kismode)
	{
		_ref_table = (double *) calloc(_nbDiffKmers,sizeof(double));
		parse_ref();
	}
	
	count_t * count_table;
 
	if(!_kismode)
	{
		count_table = (count_t *) malloc( _nbBanks * _nbDiffKmers * sizeof(count_t));
		memset(count_table,0, _nbBanks * _nbDiffKmers * sizeof(count_t));
	}
	//printf("nb diff kmers  %llu  nb banks %lu  memory allocated %llu MB\n",_nbDiffKmers,itBanks.size(), _nbBanks * _nbDiffKmers * sizeof(count_t) /1024LL/1024LL );
	
	count_t * count_table_current_bank;
	fprintf(stderr,"Launching kmer counting with offset %i step %i \n",_offset,_step);
	
	
	_progress->init ();
	
	//loop over banks
	for (size_t ii=0; ii<itBanks.size(); ii++)
	{
		
		std::vector<std::string> tabnames;
		
		if(itBanks.size() > 1)
		{
			fprintf(_outfile,"-----Bank %zu ------\n",ii);
		}
		
		Iterator<Sequence>* itSeq = itBanks[ii];
		
		//just to count nb seq :
		_nbSeq = 0;
		for (itSeq->first(); !itSeq->isDone(); itSeq->next())
		{
			Sequence& seq = itSeq->item();
			tabnames.push_back(seq.getComment());
			_nbSeq ++;
		}
		
		_resu_table = (double *) calloc(_nbSeq,sizeof(double));
		
		//ProgressIterator<Sequence>* itSeq = itBanks[ii];
		//printf("--- Counting  Bank %zu  ----\n",ii);
		
		count_table_current_bank =  & (count_table[ _nbDiffKmers * ii]);
		//loop over sequences.
		setDispatcher (  new Dispatcher (_nb_cores) );
		
		if(_dontreverse)
			getDispatcher()->iterate (itSeq,  KmerCounter<ModelDirect>(&modeld,count_table_current_bank,_nbDiffKmers,_progress,_offset,_step,_ref_table,_resu_table,_kismode), 1000);
		else
			getDispatcher()->iterate (itSeq,  KmerCounter<ModelCanonical>(&model,count_table_current_bank,_nbDiffKmers,_progress,_offset,_step,_ref_table,_resu_table,_kismode), 1000);
		
		//output result
		
		if(_kismode)
		{
			for (u_int64_t nk= 0; nk<_nbSeq; nk++ )
			{
				fprintf(_outfile,"%s\t%f\n",tabnames[nk].c_str(),_resu_table[nk]);
				
			}
		}
		
		free(_resu_table);
	}
	
	
	
	_progress->finish ();
	
	
	
	if(!_kismode)
	{
		for(uint64_t km = 0; km <_nbDiffKmers; km++ )
		{
            kmer_type km2; km2.setVal(km);
			fprintf(_outfile,"%s",model.toString(km2).c_str());
			for (size_t ii=0; ii<itBanks.size(); ii++)
			{
				fprintf(_outfile,"\t%i", count_table[IX(km,ii)]);
				
			}
			fprintf(_outfile,"\n");
			
		}
		
	}
	
	if(getParser()->saw(STR_URI_OUTPUT))
		fclose(_outfile);
	
	if(_kismode)
		free(_ref_table);
	
	if(!_kismode)
		free(count_table);
	
}








