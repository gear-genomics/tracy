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
#include "consensus.h"

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
  std::cout << "    decompose    variant calling and indel decomposition" << std::endl;
  std::cout << "    consensus    consensus for a pair of trace files" << std::endl;
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
    std::cout << " using Boost: v" << BOOST_VERSION / 100000 << "." << BOOST_VERSION / 100 % 1000 << "." << BOOST_VERSION % 100 << std::endl;
    std::cout << " using HTSlib: v" << hts_version() << std::endl;
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
    bsd();
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
  else if ((std::string(argv[1]) == "consensus")) {
    return consensus(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "assemble")) {
    return assemble(argc-1,argv+1);
  } else {
    std::cerr << "Unrecognized command " << std::string(argv[1]) << std::endl;
    return 1;
  }
}

