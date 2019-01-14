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

#ifndef SAGE_H
#define SAGE_H

#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>

#include <sdsl/suffix_arrays.hpp>

#include <htslib/faidx.h>

#include "abif.h"
#include "align.h"
#include "gotoh.h"
#include "fmindex.h"
#include "json.h"
#include "profile.h"

using namespace sdsl;

namespace tracy {

  struct SageConfig {
    uint16_t linelimit;
    uint16_t trimLeft;
    uint16_t trimRight;
    uint16_t kmer;
    uint16_t maxindel;
    uint16_t minKmerSupport;
    float pratio;
    std::string format;
    boost::filesystem::path align;
    boost::filesystem::path outfile;
    boost::filesystem::path ab;
    boost::filesystem::path genome;
  };
  
  int sage(int argc, char** argv) {

#ifdef PROFILE
    ProfilerStart("tracy.prof");
#endif
    
    SageConfig c;
  
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,g", boost::program_options::value<boost::filesystem::path>(&c.genome), "(gzipped) fasta or wildtype ab1 file")
      ("pratio,p", boost::program_options::value<float>(&c.pratio)->default_value(0.33), "peak ratio to call base")
      ("kmer,k", boost::program_options::value<uint16_t>(&c.kmer)->default_value(15), "kmer size to anchor trace")
      ("support,s",  boost::program_options::value<uint16_t>(&c.minKmerSupport)->default_value(3), "min. kmer support")
      ("maxindel,m", boost::program_options::value<uint16_t>(&c.maxindel)->default_value(1000), "max. indel size in Sanger trace")
      ("trimLeft,l", boost::program_options::value<uint16_t>(&c.trimLeft)->default_value(50), "trim size left")
      ("trimRight,r", boost::program_options::value<uint16_t>(&c.trimRight)->default_value(50), "trim size right")
      ;
    
    boost::program_options::options_description otp("Output options");
    otp.add_options()
      ("linelimit,n", boost::program_options::value<uint16_t>(&c.linelimit)->default_value(60), "alignment line length")
      ("format,f", boost::program_options::value<std::string>(&c.format)->default_value("json"), "output format [json|align]")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.json"), "output file")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
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
      std::cout << "Usage: tracy " << argv[0] << " [OPTIONS] -g genome.fa trace.ab1" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }
    if (c.maxindel < 1) c.maxindel = 1;
    
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
    int32_t ft = traceFormat(c.ab.string());
    if (ft == 0) {
      if (!readab(c.ab.string(), tr)) return -1;
    } else if (ft == 1) {
      if (!readscf(c.ab.string(), tr)) return -1;
    } else {
      std::cerr << "Unknown trace file type!" << std::endl;
      return -1;
    }
    
    // Call bases
    BaseCalls bc;
    basecall(tr, bc, c.pratio);

    // Load reference
    csa_wt<> fm_index;
    ReferenceSlice rs;
    if (!loadFMIdx(c, rs, fm_index)) return -1;
    
    // Find reference match
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Find reference match" << std::endl;
    if (!getReferenceSlice(c, fm_index, bc, rs)) return -1;

    // Create trimmed trace and reference profile
    typedef boost::multi_array<float, 2> TProfile;
    TProfile ptrace;
    createProfile(tr, bc, ptrace, c.trimLeft, c.trimRight);
    TProfile prefslice;
    createProfile(c, rs, prefslice);
    
    // Debug Profile
    //for(uint32_t i = 0; i<ptrace.shape()[0]; ++i) {
    //for(uint32_t j = 0; j<ptrace.shape()[1]; ++j) std::cerr << ptrace[i][j];
    //std::cerr << std::endl;
    //}
    
    // Semi-global alignment
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Alignment" << std::endl;
    typedef boost::multi_array<char, 2> TAlign;
    TAlign align;
    AlignConfig<true, false> semiglobal;
    DnaScore<int> sc(5, -4, -10, -1);
    gotoh(ptrace, prefslice, align, semiglobal, sc);
    
    // Debug Alignment
    //for(uint32_t i = 0; i<align.shape()[0]; ++i) {
    //for(uint32_t j = 0; j<align.shape()[1]; ++j) std::cerr << align[i][j];
    //std::cerr << std::endl;
    //}
    
    // Trim initial reference slice and extend to full trace
    trimReferenceSlice(c, align, rs);
    TProfile ptr;
    createProfile(tr, bc, ptr);
    TProfile prs;
    createProfile(c, rs, prs);
    
    // Global alignment to trimmed reference
    typedef boost::multi_array<char, 2> TAlign;
    TAlign final;
    gotoh(ptr, prs, final, semiglobal, sc);
    // Debug Alignment
    //for(uint32_t i = 0; i<final.shape()[0]; ++i) {
    //for(uint32_t j = 0; j<final.shape()[1]; ++j) std::cerr << final[i][j];
    //std::cerr << std::endl;
    //}    
    
    // Pad the trace according to alignment
    BaseCalls nbc;
    Trace ntr;
    alignmentTracePadding(final, tr, bc, ntr, nbc);
    
    // Output
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Output" << std::endl;
    if (c.format == "align") {
      // For the command line we want to show the trimmed alignment (as in Indigo)
      TAlign final1;
      std::string pritrim = trimmedSeq(bc.primary, c.trimLeft, c.trimRight);
      int32_t alTrimmedScore = gotoh(pritrim, rs.refslice, final1, semiglobal, sc);
      plotAlignment(c, final1, rs, alTrimmedScore);
    } else traceAlignJsonOut(c.outfile.string(), nbc, ntr, rs, final);


#ifdef PROFILE
    ProfilerStop();
#endif
    
    // Done
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Done." << std::endl;

    return 0;
  }

}

#endif


