/*
============================================================================
Tracy: Trace File Handling
============================================================================
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
============================================================================
Contact: Tobias Rausch (rausch@embl.de)
============================================================================
*/

#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>

#define BOOST_DISABLE_ASSERTS

#ifdef PROFILE
#include "gperftools/profiler.h"
#endif

#include "version.h"
#include "index.h"
#include "teal.h"
#include "sage.h"
#include "indigo.h"
#include "assemble.h"

using namespace tracy;


inline void
displayUsage() {
  std::cout << "Usage: tracy <command> <arguments>" << std::endl;
  std::cout << std::endl;
  std::cout << "Commands:" << std::endl;
  std::cout << std::endl;
  std::cout << "    index        index FASTA reference file" << std::endl;
  std::cout << "    basecall     basecall Chromatogram trace file" << std::endl;
  std::cout << "    align        alignment of a trace file to a genome" << std::endl;
  std::cout << "    decompose    separate a mutated and wildtype allele" << std::endl;
  std::cout << "    assemble     assemble a set of trace files" << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
}

int main(int argc, char **argv) {
  if (argc < 2) { 
    printTitle("Tracy");
    displayUsage();
    return 0;
  }
  
  if ((std::string(argv[1]) == "version") || (std::string(argv[1]) == "--version") || (std::string(argv[1]) == "--version-only") || (std::string(argv[1]) == "-v")) {
    std::cout << "Tracy version: v" << tracyVersionNumber << std::endl;
    return 0;
  }
  else if ((std::string(argv[1]) == "help") || (std::string(argv[1]) == "--help") || (std::string(argv[1]) == "-h") || (std::string(argv[1]) == "-?")) {
    printTitle("Tracy");
    displayUsage();
    return 0;
  }
  else if ((std::string(argv[1]) == "warranty") || (std::string(argv[1]) == "--warranty") || (std::string(argv[1]) == "-w")) {
    displayWarranty();
    return 0;
  }
  else if ((std::string(argv[1]) == "license") || (std::string(argv[1]) == "--license") || (std::string(argv[1]) == "-l")) {
    gplV3();
    return 0;
  }
  else if ((std::string(argv[1]) == "index")) {
    return index(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "basecall")) {
    return teal(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "align")) {
    return sage(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "decompose")) {
    return indigo(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "assemble")) {
    return assemble(argc-1,argv+1);
  } else {
    std::cerr << "Unrecognized command " << std::string(argv[1]) << std::endl;
    return 1;
  }
}

