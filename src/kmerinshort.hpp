/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#ifndef _TOOL_kis_HPP_
#define _TOOL_kis_HPP_

/********************************************************************************/
#include <gatb/gatb_core.hpp>
/********************************************************************************/


using namespace gatb::core;
using namespace gatb::core::tools;
using namespace gatb::core::bank;
using namespace gatb::core::kmer::impl;

using namespace gatb::core::system;
using namespace gatb::core::system::impl;

typedef Kmer<>::Type  kmer_type;
typedef Kmer<>::ModelCanonical    ModelCanonical;
typedef Kmer<>::ModelDirect    ModelDirect;

//typedef gatb::core::kmer::impl::Kmer<>::ModelDirect::Iterator KmerIteratorDirect;


class kis : public Tool
{
	
	static const char* STR_REFK;
	
	void parse_ref ();
	
	
public:
	bool _dontreverse;
	size_t          _kmerSize;
	size_t          _nbBanks;
	uint64_t _nbDiffKmers;
	uint64_t _nbSeq;
	int _offset;
	int _step;
	bool _kismode;
	
	IBank* _inputBank;
	std::string _inputFilename;
	std::string _kmerValFilename;
	std::string _outputFilename;
	FILE * _outfile;
	FILE * _kmerRefFile;
	
	double * _ref_table;
	double * _resu_table;
	
	int _nb_cores;
	
	gatb::core::tools::dp::IteratorListener* _progress;
	void setProgress (gatb::core::tools::dp::IteratorListener* progress)  { SP_SETATTR(progress); }
	
	
	
	
	// Constructor
	kis ();
	
	void execute ();
	
};


/********************************************************************************/

#endif

