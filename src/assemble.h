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
    float matchFraction;
    float fractionCalled;
    std::string outprefix;
    DnaScore<int32_t> aliscore;
    boost::filesystem::path reference;
    boost::filesystem::path alignment;
    std::vector<boost::filesystem::path> ab;
  };


  struct TraceScore {
    int32_t score;
    int32_t idx;
    bool forward;

    TraceScore(int32_t const s, int32_t const i, bool const f) : score(s), idx(i), forward(f) {}
  };


  template<typename TTraceScore>
  struct SortTraceScore : public std::binary_function<TTraceScore, TTraceScore, bool>
  {
    inline bool operator()(TTraceScore const& ts1, TTraceScore const& ts2) {
      return ((ts1.score > ts2.score) || ((ts1.score == ts2.score) && (ts1.idx < ts2.idx)));
    }
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
    //std::cerr << "Sequence alignment" << std::endl;
    for(TAIndex i = 0; i < (TAIndex) in.shape()[0]; ++i) {
      for(TAIndex j = 0; j < (TAIndex) in.shape()[1]; ++j) {
	out[i][j] = in[i][j];
	//std::cerr << out[i][j];
      }
      //std::cerr << std::endl;
    }
    //std::cerr << std::endl;
  }
  
  int assemble(int argc, char** argv) {
    AssembleConfig c;
  
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("pratio,p", boost::program_options::value<float>(&c.pratio)->default_value(0.33), "peak ratio to call base")
      ("trim,t", boost::program_options::value<float>(&c.trimStringency)->default_value(4), "trimming stringency [1:9]")
      ("fracmatch,f", boost::program_options::value<float>(&c.matchFraction)->default_value(0.5), "min. fraction of matches [0:1]")
      ("called,d", boost::program_options::value<float>(&c.fractionCalled)->default_value(0.1), "fraction of traces required for consensus")
      ("reference,r", boost::program_options::value<boost::filesystem::path>(&c.reference), "reference-guided assembly (optional)")
      ("outprefix,o", boost::program_options::value<std::string>(&c.outprefix)->default_value("out"), "output prefix")
      ;

    boost::program_options::options_description alignment("Alignment scoring options");
    alignment.add_options()
      ("gapopen,g", boost::program_options::value<int32_t>(&c.gapopen)->default_value(-10), "gap open")
      ("gapext,e", boost::program_options::value<int32_t>(&c.gapext)->default_value(-4), "gap extension")
      ("match,m", boost::program_options::value<int32_t>(&c.match)->default_value(3), "match")
      ("mismatch,n", boost::program_options::value<int32_t>(&c.mismatch)->default_value(-5), "mismatch")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.ab), "ab1")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(alignment).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(alignment);
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

    // Check match fraction
    if (c.matchFraction < 0) c.matchFraction = 0;
    else if (c.matchFraction > 1) c.matchFraction = 1;

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
      _createProfile(seq, prefslice);
      
      // Traces and alignment score objects
      std::vector<TProfile> traceProfiles;
      std::vector<TraceScore> scoreIdx;
      
      // Load *.ab1 files
      now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Align ab1 files" << std::endl;
      for(uint32_t i = 0; i < c.ab.size(); ++i) {
	std::cout << "Processing " << c.ab[i].string() << " [" << i << "]" << std::endl;
	
	Trace tr;
	int32_t ft = traceFormat(c.ab[i].string());
	if (ft == 0) {
	  if (!readab(c.ab[i].string(), tr)) return -1;
	} else if (ft == 1) {
	  if (!readscf(c.ab[i].string(), tr)) return -1;
	} else {
	  std::cerr << "Unknown trace file type!" << std::endl;
	  return -1;
	}

	// Call bases
	BaseCalls bc;
	basecall(tr, bc, c.pratio);

	// Get trim sizes
	uint32_t trimLeft = 0;
	uint32_t trimRight = 0;
	trimTrace(c, bc, trimLeft, trimRight);
	
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

	// Final score
	double seqsize = ptrace.shape()[1];
	double scoreThreshold = seqsize * c.matchFraction * c.aliscore.match + seqsize * (1 - c.matchFraction) * c.aliscore.mismatch; // 60% matches
	//std::cerr << scoreThreshold << ',' << gsFwd << ',' << gsRev << std::endl;
	if ((gsFwd > scoreThreshold) || (gsRev > scoreThreshold)) {
	  int32_t bestScore = std::max(gsFwd, gsRev);
	  if (gsFwd >= gsRev) {
	    scoreIdx.push_back(TraceScore(bestScore, i, true));
	    traceProfiles.push_back(ptrace);
	  } else {
	    scoreIdx.push_back(TraceScore(bestScore, i, false));
	    traceProfiles.push_back(prevtrace);
	  }
	} else {
	  std::cerr << "Warning: " << c.ab[i].string() << " is not matching to the reference! Trace file will be excluded!" << std::endl;
	  // Push-back empty trace and sequence to keep the input file order
	  TProfile empty;
	  traceProfiles.push_back(empty);
	}
      }

      // Sort score
      std::sort(scoreIdx.begin(), scoreIdx.end(), SortTraceScore<TraceScore>());

      // Align iteratively
      if (scoreIdx.size()) {
	AlignConfig<true, false> semiglobal;
	gotoh(traceProfiles[scoreIdx[0].idx], prefslice, align, semiglobal, c.aliscore);
	for(uint32_t i = 1; i < scoreIdx.size(); ++i) {
	  TAlign alignNew;
	  TProfile alignProfile;
	  _createProfile(align, alignProfile);
	  gotoh(traceProfiles[scoreIdx[i].idx], alignProfile, alignNew, semiglobal, c.aliscore);

	  // Debug profile alignment
	  //std::cerr << "Profile alignment" << std::endl;
	  //for(uint32_t i = 0; i < alignNew.shape()[0]; ++i) {
	  //for(uint32_t j = 0; j < alignNew.shape()[1]; ++j) std::cerr << alignNew[i][j];
	  //std::cerr << std::endl;
	  //}
	  
	  // Create new sequence alignment based on profile alignment
	  TAlign alignCombined;
	  uint32_t nSeq = align.shape()[0] + 1;
	  uint32_t nCol = alignNew.shape()[1];
	  uint32_t ap = 0;
	  alignCombined.resize(boost::extents[nSeq][alignNew.shape()[1]]);
	  for(uint32_t j = 0; j < nCol; ++j) {
	    alignCombined[0][j] = alignNew[0][j];
	    if (alignNew[1][j] != '-') {
	      for(uint32_t k = 1; k < nSeq; ++k) alignCombined[k][j] = align[k-1][ap];
	      ++ap;
	    } else {
	      for(uint32_t k = 1; k < nSeq; ++k) alignCombined[k][j] = '-';
	    }
	  }

	  // Overwrite old alignment
	  overwriteArray(alignCombined, align);
	}

	// Consensus calling
	consensus(c, align, gapped, cs);

	// Output horizontal alignment
	std::string alignfilename = c.outprefix + ".align.fa";
	std::ofstream vfile(alignfilename.c_str());
	typedef typename TAlign::index TAIndex;
	for(uint32_t i = 0; i < scoreIdx.size(); ++i) {
	  int32_t alignRow = scoreIdx.size()-i-1;
	  vfile	<< ">" << c.ab[scoreIdx[i].idx].stem().string() << std::endl;
	  for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
	    vfile << align[alignRow][j];
	  }
	  vfile << std::endl;
	}
	vfile << ">Reference" << std::endl;
	for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
	  vfile << align[scoreIdx.size()][j];
	}
	vfile << std::endl;
	vfile.close();
	
	// Output MSA and gapped traces
	std::string filename = c.outprefix + ".json";
	std::ofstream rfile(filename.c_str());
	rfile << "{" << std::endl;
	rfile << "\"gapFreeConsensus\": \"" << cs << "\"," << std::endl;
	rfile << "\"gappedConsensus\": \"" << gapped << "\"," << std::endl;
	rfile << "\"msa\": " << std::endl;
	rfile << "[" << std::endl;
	for(uint32_t i = 0; i < scoreIdx.size(); ++i) {
	  if (i!=0) rfile << ',' << std::endl;
	  alignedTraceByRow(rfile, align, scoreIdx.size()-i-1, c.ab[scoreIdx[i].idx].stem().string(), scoreIdx[i].forward, false);
	}
	rfile << ',' << std::endl;
	alignedTraceByRow(rfile, align, scoreIdx.size(), "", true, true);
	rfile << "]," << std::endl;
	rfile << "\"gappedTraces\": " << std::endl;
	rfile << "[" << std::endl;
	for(uint32_t i = 0; i < scoreIdx.size(); ++i) {
	  if (i!=0) rfile << ", ";
	  Trace tr;
	  int32_t ft = traceFormat(c.ab[i].string());
	  if (ft == 0) {
	    if (!readab(c.ab[i].string(), tr)) return -1;
	  } else if (ft == 1) {
	    if (!readscf(c.ab[i].string(), tr)) return -1;
	  } else {
	    std::cerr << "Unknown trace file type!" << std::endl;
	    return -1;
	  }

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

	  // Reverse complement trace if necessary
	  BaseCalls padbc;
	  Trace padtr;
	  if (scoreIdx[i].forward) {
	    // Debug
	    //std::string filename = c.ab[scoreIdx[i].idx].stem().string() + ".txt";
	    //traceTxtOut(filename, nbc, ntr);
	    alignmentTracePadding(align, ntr, nbc, scoreIdx.size()-i-1, padtr, padbc);
	  } else {
	    BaseCalls tbc;
	    Trace ttr;
	    // Reverese complement
	    reverseComplementTrace(ntr, nbc, ttr, tbc);
	    alignmentTracePadding(align, ttr, tbc, scoreIdx.size()-i-1, padtr, padbc);
	  }

	  // Append gapped trace to output
	  assemblyTrace(rfile, padbc, padtr, c.ab[scoreIdx[i].idx].stem().string());
	}
	rfile << "]" << std::endl;
	rfile << "}" << std::endl;
	rfile.close();
      }
    } else {
      // De-novo assembly
      std::cout << "Please specify a reference file!" << std::endl;
      return -1;
      
      // Load *.ab1 files
      now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Load ab1 files" << std::endl;
      std::vector<SequenceSegment> seqSegment;
      for(uint32_t i = 0; i < c.ab.size(); ++i) {
	Trace tr;
	int32_t ft = traceFormat(c.ab[i].string());
	if (ft == 0) {
	  if (!readab(c.ab[i].string(), tr)) return -1;
	} else if (ft == 1) {
	  if (!readscf(c.ab[i].string(), tr)) return -1;
	} else {
	  std::cerr << "Unknown trace file type!" << std::endl;
	  return -1;
	}

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
    std::string filename = c.outprefix + ".vertical";
    std::ofstream vfile(filename.c_str());
    typedef typename TAlign::index TAIndex;
    for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
      for(TAIndex i = 0; i < (TAIndex) align.shape()[0]; ++i) {
	vfile << align[i][j];
      }
      vfile << '|' << gapped[j] << std::endl;
    }
    vfile.close();

    // Done
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Done." << std::endl;
    return 0;
  }

}

#endif
