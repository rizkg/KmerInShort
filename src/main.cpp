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

// We include the header file for the tool
#include <kmerinshort.hpp>

/********************************************************************************/

int main (int argc, char* argv[])
{
	
	if(argc==2 &&  (strcmp(argv[1],"--version")==0  ||  strcmp(argv[1],"-v")==0) ){
		printf("KmerInShort version 1.1.1\n");
		return EXIT_SUCCESS;
	}
	
    try
    {
        // We run the tool with the provided command line arguments.
        kis().run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

