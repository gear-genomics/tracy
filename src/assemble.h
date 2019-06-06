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

#ifndef ASSEMBLE_H
#define ASSEMBLE_H

#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#include "msa.h"
#include "trim.h"

using namespace sdsl;

namespace tracy {
  
  struct AssembleConfig {
    bool hasReference;
    int32_t gapopen;
    int32_t gapext;
    int32_t match;
    int32_t mismatch;
    float pratio;
    float trimStringency;
    float fractionCalled;
    DnaScore<int32_t> aliscore;
    boost::filesystem::path reference;
    boost::filesystem::path alignment;
    boost::filesystem::path outfile;
    std::vector<boost::filesystem::path> ab;
  };

  struct SequenceSegment {
    std::string seq;
    int32_t trimLeft;
    int32_t trimRight;
    bool forward;
    SequenceSegment(std::string const& s, int32_t tl, int32_t tr, bool f) : seq(s), trimLeft(tl), trimRight(tr), forward(f) {}
  };


  template<typename TArray>
  inline void
  overwriteArray(TArray const& in, TArray& out) {
    typedef typename TArray::index TAIndex;
    out.resize(boost::extents[in.shape()[0]][in.shape()[1]]);
    for(TAIndex i = 0; i < (TAIndex) in.shape()[0]; ++i) {
      for(TAIndex j = 0; j < (TAIndex) in.shape()[1]; ++j) {
	out[i][j] = in[i][j];
      }
    }
  }
  
  int assemble(int argc, char** argv) {
    AssembleConfig c;
  
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("pratio,p", boost::program_options::value<float>(&c.pratio)->default_value(0.33), "peak ratio to call base")
      ("trim,t", boost::program_options::value<float>(&c.trimStringency)->default_value(4), "trimming stringency [1:9]")
      ("called,d", boost::program_options::value<float>(&c.fractionCalled)->default_value(0.1), "fraction of traces required for consensus")
      ("reference,r", boost::program_options::value<boost::filesystem::path>(&c.reference), "reference-guided assembly (optional)")
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

    // Check reference
    if (vm.count("reference")) c.hasReference = true;
    else c.hasReference = false;

    // Alignment Scoring
    c.aliscore = DnaScore<int32_t>(c.match, c.mismatch, c.gapopen, c.gapext);
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "tracy ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;

    // Multiple sequence alignment
    typedef boost::multi_array<char, 2> TAlign;
    TAlign align;
    std::string gapped;
    std::string cs;
    
    // Reference-guided
    if (c.hasReference) {
      // Load reference
      std::string faname = "";
      std::string seq = "";
      loadSingleFasta(c.reference.string(), faname, seq);

      // Reference profile
      typedef boost::multi_array<float, 2> TProfile;
      TProfile prefslice;
      createProfile(seq, prefslice);
      
      // Traces and alignment score objects
      std::vector<TProfile> traceProfiles;
      typedef std::pair<int32_t, int32_t> TScoreIdx;
      std::vector<TScoreIdx> scoreIdx;
      std::vector<std::string> sequences;
      
      // Load *.ab1 files
      now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Align ab1 files" << std::endl;
      for(uint32_t i = 0; i < c.ab.size(); ++i) {
	Trace tr;
	if (!readab(c.ab[i].string(), tr)) return -1;

	// Call bases
	BaseCalls bc;
	basecall(tr, bc, c.pratio);

	// Get trim sizes
	uint32_t trimLeft = 0;
	uint32_t trimRight = 0;
	trimTrace(c, bc, trimLeft, trimRight);
	std::string primarySeq = trimmedSeq(bc.primary, trimLeft, trimRight);
	
	// Create Trace Profile
	TProfile ptrace;
	createProfile(tr, bc, ptrace, trimLeft, trimRight);
	
	// Align
	AlignConfig<true, false> semiglobal;
	int32_t gsFwd = gotohScore(ptrace, prefslice, semiglobal, c.aliscore);

	// Reverse complement profile
	TProfile prevtrace;
	reverseComplementProfile(ptrace, prevtrace);
	int32_t gsRev = gotohScore(prevtrace, prefslice, semiglobal, c.aliscore);

	// Debug
	std::cerr << gsFwd << ',' << gsRev << ',' << c.ab[i] << std::endl;
	
	// Final score
	double seqsize = ptrace.shape()[1];
	double scoreThreshold = seqsize * 0.6 * c.aliscore.match + seqsize * 0.4 * c.aliscore.mismatch; // 60% matches
	if ((gsFwd > scoreThreshold) || (gsRev > scoreThreshold)) {
	  int32_t bestScore = std::max(gsFwd, gsRev);
	  scoreIdx.push_back(std::make_pair(-bestScore, i));
	  if (gsFwd >= gsRev) {
	    std::cerr << "Forward alignment" << std::endl;
	    traceProfiles.push_back(ptrace);
	    sequences.push_back(primarySeq);
	  } else {
	    std::cerr << "Reverse alignment" << std::endl;
	    traceProfiles.push_back(prevtrace);
	    reverseComplement(primarySeq);
	    sequences.push_back(primarySeq);
	  }
	} else {
	  std::cerr << "Warning: " << c.ab[i].string() << " is not matching to the reference! Trace file will be excluded!" << std::endl;
	  // Push-back empty trace to keep the input file order
	  TProfile empty;
	  traceProfiles.push_back(empty);
	}
      }

      // Sort score
      std::sort(scoreIdx.begin(), scoreIdx.end());

      // Align iteratively
      if (scoreIdx.size()) {
	AlignConfig<true, false> semiglobal;
	gotoh(traceProfiles[scoreIdx[0].second], prefslice, align, semiglobal, c.aliscore);
	for(uint32_t i = 1; i < scoreIdx.size(); ++i) {
	  TAlign alignSeq;
	  alignSeq.resize(boost::extents[1][sequences[scoreIdx[i].second].size()]);
	  uint32_t ind = 0;
	  for(typename std::string::const_iterator str = sequences[scoreIdx[i].second].begin(); str != sequences[scoreIdx[i].second].end(); ++str) alignSeq[0][ind++] = *str;
	  TAlign alignNew;	  
	  gotoh(alignSeq, align, alignNew, semiglobal, c.aliscore);
	  overwriteArray(alignNew, align);
	}

	// Consensus calling
	consensus(c, align, gapped, cs);

	// Output MSA and gapped traces
	std::ofstream rfile(c.outfile.c_str());
	rfile << "{" << std::endl;
	rfile << "\"gapFreeConsensus\": \"" << cs << "\"," << std::endl;
	rfile << "\"gappedConsensus\": \"" << gapped << "\"," << std::endl;
	rfile << "\"msa\": " << std::endl;
	rfile << "[" << std::endl;
	for(uint32_t i = 0; i < scoreIdx.size(); ++i) {
	  if (i!=0) rfile << ',' << std::endl;
	  alignedTraceByRow(rfile, align, scoreIdx.size()-i-1, c.ab[scoreIdx[i].second].stem().string(), false);
	}
	rfile << ',' << std::endl;
	alignedTraceByRow(rfile, align, scoreIdx.size(), "", true);
	rfile << "]," << std::endl;
	rfile << "\"gappedTraces\": " << std::endl;
	rfile << "[" << std::endl;
	for(uint32_t i = 0; i < scoreIdx.size(); ++i) {
	  if (i!=0) rfile << ", ";
	  Trace tr;
	  if (!readab(c.ab[scoreIdx[i].second].string(), tr)) return -1;

	  // Call bases
	  BaseCalls bc;
	  basecall(tr, bc, c.pratio);

	  // Get trim sizes
	  uint32_t trimLeft = 0;
	  uint32_t trimRight = 0;
	  trimTrace(c, bc, trimLeft, trimRight);

	  // Hard trim of the trace data structures
	  BaseCalls nbc;
	  Trace ntr;
	  trimTrace(tr, bc, trimLeft, trimRight, ntr, nbc);

	  // Debug
	  //std::string filename = c.ab[scoreIdx[i].second].stem().string() + ".txt";
	  //traceTxtOut(filename, nbc, ntr);

	  // Trace padding with gaps
	  BaseCalls padbc;
	  Trace padtr;
	  alignmentTracePadding(align, ntr, nbc, scoreIdx.size()-i-1, padtr, padbc);

	  // Append gapped trace to output
	  assemblyTrace(rfile, padbc, padtr, c.ab[scoreIdx[i].second].stem().string());
	}
	rfile << "]" << std::endl;
	rfile << "}" << std::endl;
	rfile.close();
      }
    } else {
      // De-novo assembly
      
      // Load *.ab1 files
      now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Load ab1 files" << std::endl;
      std::vector<SequenceSegment> seqSegment;
      for(uint32_t i = 0; i < c.ab.size(); ++i) {
	Trace tr;
	if (!readab(c.ab[i].string(), tr)) return -1;

	// Call bases
	BaseCalls bc;
	basecall(tr, bc, c.pratio);

	// Get trim sizes
	uint32_t trimLeft = 0;
	uint32_t trimRight = 0;
	trimTrace(c, bc, trimLeft, trimRight);
	seqSegment.push_back(SequenceSegment(bc.primary, trimLeft, trimRight, true));
      }

      // Optimize layout/trimming
      now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Optimize layout/trimming" << std::endl;
      std::vector<std::string> traceSet;
      revSeqBasedOnDist(c, seqSegment, traceSet);	

      // Assemble
      now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Assemble traces" << std::endl;
      std::string consensus;
      msa(c, traceSet, consensus);
    }

    // Output vertical alignment
    boost::iostreams::filtering_ostream rcfile;
    rcfile.push(boost::iostreams::gzip_compressor());
    rcfile.push(boost::iostreams::file_sink(c.alignment.c_str(), std::ios_base::out | std::ios_base::binary));
    typedef typename TAlign::index TAIndex;
    for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
      for(TAIndex i = 0; i < (TAIndex) align.shape()[0]; ++i) {
	rcfile << align[i][j];
      }
      rcfile << '|' << gapped[j] << std::endl;
    }
    rcfile.pop();
    
    // Done
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Done." << std::endl;
    return 0;
  }

}

#endif
