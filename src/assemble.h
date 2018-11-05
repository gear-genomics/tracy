/*
============================================================================
Tracy: Trace File Handling
============================================================================
Copyright (C) 2017,2018 Tobias Rausch

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

#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#include "msa.h"
#include "trim.h"

using namespace sdsl;

namespace tracy {
  
  struct AssembleConfig {
    uint16_t linelimit;
    int32_t gapopen;
    int32_t gapext;
    int32_t match;
    int32_t mismatch;
    float pratio;
    float trimStringency;
    float fractionCalled;
    std::string format;
    DnaScore<int32_t> aliscore;
    boost::filesystem::path alignment;
    boost::filesystem::path outfile;
    std::vector<boost::filesystem::path> ab;
  };

  int assemble(int argc, char** argv) {
    AssembleConfig c;
  
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("pratio,p", boost::program_options::value<float>(&c.pratio)->default_value(0.33), "peak ratio to call base")
      ("trim,t", boost::program_options::value<float>(&c.trimStringency)->default_value(4), "trimming stringency [1:9]")
      ("called,d", boost::program_options::value<float>(&c.fractionCalled)->default_value(0.1), "fraction of traces required for consensus")
      ;

    boost::program_options::options_description alignment("Alignment scoring options");
    alignment.add_options()
      ("gapopen,g", boost::program_options::value<int32_t>(&c.gapopen)->default_value(-10), "gap open")
      ("gapext,e", boost::program_options::value<int32_t>(&c.gapext)->default_value(-4), "gap extension")
      ("match,m", boost::program_options::value<int32_t>(&c.match)->default_value(3), "match")
      ("mismatch,n", boost::program_options::value<int32_t>(&c.mismatch)->default_value(-5), "mismatch")
      ;
    
    boost::program_options::options_description otp("Output options");
    otp.add_options()
      ("linelimit,l", boost::program_options::value<uint16_t>(&c.linelimit)->default_value(60), "alignment line length")
      ("format,f", boost::program_options::value<std::string>(&c.format)->default_value("json"), "output format [json|align]")
      ("alignment,a", boost::program_options::value<boost::filesystem::path>(&c.alignment)->default_value("al.fa.gz"), "vertical alignment")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.json"), "output file")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.ab), "ab1")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(alignment).add(otp).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(alignment).add(otp);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    
    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file"))) {
      std::cout << "Usage: tracy " << argv[0] << " [OPTIONS] trace1.ab1 trace2.ab1 ..." << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    // Check ab1
    for(uint32_t i = 0; i < c.ab.size(); ++i) {
      if (!(boost::filesystem::exists(c.ab[i]) && boost::filesystem::is_regular_file(c.ab[i]) && boost::filesystem::file_size(c.ab[i]))) {
	std::cerr << "Trace file is missing: " << c.ab[i].string() << std::endl;
	return 1;
      }
    }

    // Alignment Scoring
    c.aliscore = DnaScore<int32_t>(c.match, c.mismatch, c.gapopen, c.gapext);
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "tracy ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    // Load *.ab1 files
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Load ab1 files" << std::endl;
    std::vector<std::string> traceSet;
    for(uint32_t i = 0; i < c.ab.size(); ++i) {
      Trace tr;
      if (!readab(c.ab[i].string(), tr)) return -1;

      // Call bases
      BaseCalls bc;
      basecall(tr, bc, c.pratio);

      // Append to trace set
      uint32_t trimLeft = 0;
      uint32_t trimRight = 0;
      trimTrace(c, bc, trimLeft, trimRight);
      if (trimLeft + trimRight < bc.primary.size()) {
	std::cout << c.ab[i].string() << "(Size: " << bc.primary.size() << ", Trimming Left: " << trimLeft << ", Trimming Right: " << trimRight << ")" << std::endl;
	int32_t sz = bc.primary.size() - trimLeft - trimRight;
	traceSet.push_back(bc.primary.substr(trimLeft, sz));
      }
    }

    // Assemble Traces
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Assemble traces" << std::endl;
    std::string consensus;
    msa(c, traceSet, consensus);
    
    // Done
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Done." << std::endl;
    return 0;
  }

}

#endif
