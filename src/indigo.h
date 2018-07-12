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

#ifndef INDIGO_H
#define INDIGO_H

#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#include "decompose.h"

using namespace sdsl;

namespace tracy {
  
  struct IndigoConfig {
    uint16_t linelimit;
    uint16_t trimLeft;
    uint16_t trimRight;
    uint16_t kmer;
    uint16_t maxindel;
    uint16_t madc;
    float pratio;
    std::string outprefix;
    boost::filesystem::path ab;
    boost::filesystem::path genome;
  };

  int indigo(int argc, char** argv) {
    IndigoConfig c;
  
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "(gzipped) fasta or wildtype ab1 file")
      ("pratio,p", boost::program_options::value<float>(&c.pratio)->default_value(0.33), "peak ratio to call base")
      ("kmer,k", boost::program_options::value<uint16_t>(&c.kmer)->default_value(15), "kmer size")
      ("maxindel,m", boost::program_options::value<uint16_t>(&c.maxindel)->default_value(1000), "max. indel size in Sanger trace")
      ("trimLeft,l", boost::program_options::value<uint16_t>(&c.trimLeft)->default_value(50), "trim size left")
      ("trimRight,r", boost::program_options::value<uint16_t>(&c.trimRight)->default_value(50), "trim size right")
      ;
    
    boost::program_options::options_description otp("Output options");
    otp.add_options()
      ("output,o", boost::program_options::value<std::string>(&c.outprefix)->default_value("align"), "output file prefix")
      ("linelimit,l", boost::program_options::value<uint16_t>(&c.linelimit)->default_value(60), "alignment line length")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("madc,c", boost::program_options::value<uint16_t>(&c.madc)->default_value(5), "MAD cutoff")
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.ab), "ab1")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(otp).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(otp);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    
    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file"))) {
      std::cout << "Usage: tracy " << argv[0] << " [OPTIONS] trace.ab1" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }
    if (c.maxindel < 1) c.maxindel = 1;

    // Check ab1
    if (!(boost::filesystem::exists(c.ab) && boost::filesystem::is_regular_file(c.ab) && boost::filesystem::file_size(c.ab))) {
      std::cerr << "Trace file is missing: " << c.ab.string() << std::endl;
      return 1;
    }
    
    // Check reference
    if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
      std::cerr << "Reference file is missing: " << c.genome.string() << std::endl;
      return 1;
    }
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "tracy ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    // Load *.ab1 file
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Load ab1 file" << std::endl;
    Trace tr;
    if (!readab(c.ab.string(), tr)) return -1;

    // Call bases
    BaseCalls bc;
    basecall(tr, bc, c.pratio);

    // Create trimmed trace profile
    typedef boost::multi_array<double, 2> TProfile;
    TProfile ptrace;
    createProfile(tr, bc, ptrace, c.trimLeft, c.trimRight);

    // Identify position of indel shift in Sanger trace
    TraceBreakpoint bp;
    findBreakpoint(ptrace, bp);
    //std::cerr << "Breakpoint: " << bp.indelshift << ',' << bp.traceleft << ',' << bp.breakpoint << ',' << bp.bestDiff << std::endl;
    
    // Load reference
    csa_wt<> fm_index;
    ReferenceSlice rs;
    if (!loadFMIdx(c, rs, fm_index)) return -1;
      
    // Find reference match
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Find Reference Match" << std::endl;
    
    // Get reference slice
    if (!getReferenceSlice(c, fm_index, bc, rs)) return -1;

    // Create reference profile
    TProfile prefslice;
    createProfile(c, rs, prefslice);

    // Semi-global alignment
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Alignment" << std::endl;
    typedef boost::multi_array<char, 2> TAlign;
    TAlign align;
    AlignConfig<true, false> semiglobal;
    DnaScore<int> sc(5, -4, -10, -1);
    gotoh(ptrace, prefslice, align, semiglobal, sc);
    
    // Do we have a shifted trace?
    if (!bp.indelshift) {
      // Find breakpoint for hom. indels
      now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Homozygous InDel Search" << std::endl;
      if (!findHomozygousBreakpoint(align, bp)) return -1;
    }

    // Debug Breakpoint & Alignment
    std::cerr << "Breakpoint: " << bp.indelshift << ',' << bp.traceleft << ',' << bp.breakpoint << ',' << bp.bestDiff << std::endl;
    for(uint32_t i = 0; i<align.shape()[0]; ++i) {
      uint32_t alignedNuc = 0;
      for(uint32_t j = 0; j<align.shape()[1]; ++j) {
	if (align[0][j] != '-') {
	  ++alignedNuc;
	  if (alignedNuc == bp.breakpoint) std::cerr << "#####";
	}
	std::cerr << align[i][j];
      }
      std::cerr << std::endl;
    }    

    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Decompose Chromatogram" << std::endl;
    
    // Decompose alleles
    if (!decomposeAlleles(c, align, bc, bp, rs)) return -1;

    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Align to reference" << std::endl;
    
    // Plot alignments
    typedef boost::multi_array<char, 2> TAlign;
    TAlign alignPrimary;
    DnaScore<int> scSeq(5, -4, -50, 0);
    std::string pri = trimmedSeq(bc.primary, c.trimLeft, c.trimRight);
    gotohString(pri, rs.refslice, alignPrimary, semiglobal, scSeq);
    plotAlignmentDeprecated(c, alignPrimary, rs, 1);
    
    typedef boost::multi_array<char, 2> TAlign;
    TAlign alignSecondary;
    std::string sec = trimmedSeq(bc.secondary, c.trimLeft, c.trimRight);
    gotohString(sec, rs.refslice, alignSecondary, semiglobal, scSeq);
    plotAlignmentDeprecated(c, alignSecondary, rs, 2);
      
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Done." << std::endl;
    return 0;
  }

}

#endif
