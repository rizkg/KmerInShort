[![Build Status](https://travis-ci.org/rizkg/KmerInShort.svg?branch=master)](https://travis-ci.org/rizkg/KmerInShort)

# KmerInShort
KmerInShort counts kmers from a fasta/fastq file or list of files, and outputs results in a text file. It is limited to short kmers (k<15).
It is a part of the [FEELnc](https://github.com/tderrien/FEELnc) pipeline from V.Wucher, F.Legai and T.Derrien, a pipeline to annotate long non-coding RNAs.

# Installation
To retrieve KmerInShort and its submodule (gatb-core), type 

    git clone --recursive https://github.com/rizkg/KmerInShort

Then build the tool with 

    mkdir build;  cd build;  cmake ..;  make -j 8
  

# Author
Guillaume Rizk
guillaume.rizk@algorizk.com
