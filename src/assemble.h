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
    bool incCons;
    bool incRef;
    int32_t gapopen;
    int32_t gapext;
    int32_t match;
    int32_t mismatch;
    float pratio;
    float trimStringency;
    float matchFraction;
    float fractionCalled;
    std::string outprefix;
    std::string format;
    DnaScore<int32_t> aliscore;
    boost::filesystem::path reference;
    boost::filesystem::path alignment;
    std::vector<boost::filesystem::path> ab;
  };


  struct TraceScore {
    int32_t score;
    int32_t idx;
    int32_t newidx;
    bool forward;

    TraceScore(int32_t const s, int32_t const i, int32_t const n, bool const f) : score(s), idx(i), newidx(n), forward(f) {}

    bool operator<(const TraceScore& ts2) const {
      return ((score > ts2.score) || ((score == ts2.score) && (idx < ts2.idx)));
    }
  };


  struct SequenceSegment {
    std::string seq;
    int32_t trimLeft;
    int32_t trimRight;
    bool forward;
    SequenceSegment(std::string const& s, int32_t tl, int32_t tr, bool f) : seq(s), trimLeft(tl), trimRight(tr), forward(f) {}
  };


  int assemble(int argc, char** argv) {
    AssembleConfig c;
  
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("reference,r", boost::program_options::value<boost::filesystem::path>(&c.reference), "reference-guided assembly (optional)")
      ("pratio,p", boost::program_options::value<float>(&c.pratio)->default_value(0.33), "peak ratio to call base")
      ("trim,t", boost::program_options::value<float>(&c.trimStringency)->default_value(4), "trimming stringency [1:9], 0: disable trimming")
      ("fracmatch,f", boost::program_options::value<float>(&c.matchFraction)->default_value(0.5), "min. fraction of matches [0:1]")
      ;

    boost::program_options::options_description alignment("Alignment options");
    alignment.add_options()
      ("gapopen,g", boost::program_options::value<int32_t>(&c.gapopen)->default_value(-10), "gap open")
      ("gapext,e", boost::program_options::value<int32_t>(&c.gapext)->default_value(-4), "gap extension")
      ("match,m", boost::program_options::value<int32_t>(&c.match)->default_value(3), "match")
      ("mismatch,n", boost::program_options::value<int32_t>(&c.mismatch)->default_value(-5), "mismatch")
      ;

    boost::program_options::options_description otp("Output options");
    otp.add_options()
      ("called,d", boost::program_options::value<float>(&c.fractionCalled)->default_value(0.1), "fraction of traces required for consensus")
      ("outprefix,o", boost::program_options::value<std::string>(&c.outprefix)->default_value("out"), "output prefix")
      ("format,a", boost::program_options::value<std::string>(&c.format)->default_value("fasta"), "consensus output format [fasta|fastq]")
      ("inccons,i", "include consensus in FASTA align")
      ("incref,j", "include reference in consensus computation (req. --reference)")
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

    // Check include consensus flag
    if (vm.count("inccons")) c.incCons = true;
    else c.incCons = false;

    // Check include reference flag
    if (vm.count("incref")) c.incRef = true;
    else c.incRef = false;

    // Check trimming stringency
    if (c.trimStringency != 0) {
      if (c.trimStringency > 9) c.trimStringency = 9;
      else if (c.trimStringency < 1) c.trimStringency = 1;
    }
    
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

    // Sequence Profile
    typedef boost::multi_array<float, 2> TProfile;

    // Multiple sequence alignment
    typedef boost::multi_array<char, 2> TAlign;
    TAlign align;
    std::string gapped;
    std::string cs;
    std::string qstr;
    
    // Reference-guided
    if (c.hasReference) {
      // Load reference
      std::string faname = "";
      std::string seq = "";
      if (!loadSingleFasta(c.reference.string(), faname, seq)) return -1;

      // Check reference size
      if (seq.size() > MAX_SINGLE_FASTA_SIZE) {
	std::cerr << "Reference is larger than 50Kbp. Please use a smaller reference slice!" << std::endl;
	return -1;
      }

      // Reference profile
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
	if (c.trimStringency) {
	  trimTrace(c, bc, trimLeft, trimRight);
	  if (trimLeft + trimRight >= bc.bcPos.size()) {
	    std::cerr << "Too stringent trimming parameters!" << std::endl;
	    std::cerr << c.ab[i].string() << " has no peaks left!" << std::endl;
	    return -1;
	  }
	}
	
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
	    scoreIdx.push_back(TraceScore(bestScore, i, scoreIdx.size(), true));
	    traceProfiles.push_back(ptrace);
	  } else {
	    scoreIdx.push_back(TraceScore(bestScore, i, scoreIdx.size(), false));
	    traceProfiles.push_back(prevtrace);
	  }
	} else {
	  std::cerr << "Warning: " << c.ab[i].stem().string() << " is not matching to the reference! Trace file will be excluded!" << std::endl;
	}
      }

      // Sort score
      std::sort(scoreIdx.begin(), scoreIdx.end());

      // Align iteratively
      if (scoreIdx.size()) {
	AlignConfig<true, false> semiglobal;
	gotoh(traceProfiles[scoreIdx[0].newidx], prefslice, align, semiglobal, c.aliscore);
	for(uint32_t i = 1; i < scoreIdx.size(); ++i) {
	  TAlign alignNew;
	  TProfile alignProfile;
	  _createProfile(align, alignProfile);
	  gotoh(traceProfiles[scoreIdx[i].newidx], alignProfile, alignNew, semiglobal, c.aliscore);

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
	if (c.incRef) consensus(c, align, gapped, cs, qstr, false);
	else consensus(c, align, gapped, cs, qstr, true);

	// Output horizontal alignment
	std::string alignfilename = c.outprefix + ".align.fa";
	std::ofstream vfile(alignfilename.c_str());
	typedef typename TAlign::index TAIndex;
	for(uint32_t i = 0; i < scoreIdx.size(); ++i) {
	  int32_t alignRow = scoreIdx.size()-i-1;
	  vfile	<< ">" << c.ab[scoreIdx[i].idx].stem().string();
	  if (scoreIdx[i].forward) vfile << " (forward)" << std::endl;
	  else vfile << " (reverse)" << std::endl;
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
	if (c.incCons) {
	  vfile << ">Consensus" << std::endl;
	  vfile << gapped << std::endl;
	}
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
	  int32_t ft = traceFormat(c.ab[scoreIdx[i].idx].string());
	  if (ft == 0) {
	    if (!readab(c.ab[scoreIdx[i].idx].string(), tr)) return -1;
	  } else if (ft == 1) {
	    if (!readscf(c.ab[scoreIdx[i].idx].string(), tr)) return -1;
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
	  if (c.trimStringency) trimTrace(c, bc, trimLeft, trimRight);
	  
	  // Hard trim of the basecalls, keep trace data structure
	  BaseCalls nbc;
	  trimTrace(tr, bc, trimLeft, trimRight, nbc);

	  // Reverse complement trace if necessary
	  BaseCalls padbc;
	  Trace padtr;
	  if (scoreIdx[i].forward) {
	    alignmentTracePadding(align, tr, nbc, scoreIdx.size()-i-1, padtr, padbc);
	  } else {
	    BaseCalls tbc;
	    Trace ttr;
	    // Reverese complement
	    reverseComplementTrace(tr, nbc, ttr, tbc);
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
      now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Load ab1 files" << std::endl;
      std::vector<TProfile> inputProfiles;
      std::vector<bool> fwdProfiles;
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
	if (c.trimStringency) {
	  trimTrace(c, bc, trimLeft, trimRight);
	  if (trimLeft + trimRight >= bc.bcPos.size()) {
	    std::cerr << "Too stringent trimming parameters!" << std::endl;
	    std::cerr << c.ab[i].string() << " has no peaks left!" << std::endl;
	    return -1;
	  }
	}
	
	// Create Trace Profile
	TProfile ptrace;
	createProfile(tr, bc, ptrace, trimLeft, trimRight);
	inputProfiles.push_back(ptrace);
	fwdProfiles.push_back(true);
      }

      // Optimize layout/trimming
      now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Optimize layout" << std::endl;
      revSeqBasedOnDist(c, inputProfiles, fwdProfiles);

      // Any trace we need to exclude?
      std::vector<TProfile> seqProfiles;
      std::vector<uint32_t> idxMap;
      std::vector<bool> fwd;
      for (uint32_t i = 0; i<inputProfiles.size(); ++i) {
	int32_t seqSize = inputProfiles[i].shape()[1];
	bool foundHit = false;
	for(uint32_t j = 0; j<inputProfiles.size(); ++j) {
	  if (i!=j) {
	    AlignConfig<true, true> alignconf;
	    TAlign seqAlign;
	    int32_t gs = gotoh(inputProfiles[i], inputProfiles[j], seqAlign, alignconf, c.aliscore);
	    int32_t numAligned = 0;
	    for(uint32_t k = 0; k < seqAlign.shape()[1];++k) {
	      if ((seqAlign[0][k] != '-') && (seqAlign[1][k] != '-')) ++numAligned;
	    }
	    double frac = (double) numAligned / (double) seqSize;
	    double scoreThreshold = numAligned * c.matchFraction * c.aliscore.match + numAligned * (1 - c.matchFraction) * c.aliscore.mismatch;
	    // At least 10% overlap
	    if ((frac > 0.1) && (numAligned > 25) && (gs > scoreThreshold)) {
	      foundHit = true;
	      break;
	    }
	  }
	}
	if (!foundHit) {
	  std::cerr << "Warning: " << c.ab[i].stem().string() << " is not matching to any of the other traces! Trace file will be excluded!" << std::endl;
	} else {
	  seqProfiles.push_back(inputProfiles[i]);
	  idxMap.push_back(i);
	  fwd.push_back(fwdProfiles[i]);
	}
      }

      // Enough sequences left?
      if (fwd.size() < 2) {
	std::cerr << "At least 2 traces are required for de novo assembly!" << std::endl;
	return -1;
      }

      // Assemble
      now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Assemble traces" << std::endl;
      std::vector<uint32_t> seqidx;
      msa(c, seqProfiles, align, seqidx);

      // Consensus calling
      consensus(c, align, gapped, cs, qstr, false);

      // Output horizontal alignment
      std::string alignfilename = c.outprefix + ".align.fa";
      std::ofstream vfile(alignfilename.c_str());
      typedef typename TAlign::index TAIndex;
      for(TAIndex i = 0; i < (TAIndex) align.shape()[0]; ++i) {
	vfile << ">" << c.ab[idxMap[seqidx[i]]].stem().string();
	if (fwd[seqidx[i]]) vfile << " (forward)" << std::endl;
	else vfile << " (reverse)" << std::endl;
	for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
	  vfile << align[i][j];
	}
	vfile << std::endl;
      }
      if (c.incCons) {
	vfile << ">Consensus" << std::endl;
	vfile << gapped << std::endl;
      }
      vfile.close();
      
      // Output MSA and gapped traces
      std::string filename = c.outprefix + ".json";
      std::ofstream rfile(filename.c_str());
      rfile << "{" << std::endl;
      rfile << "\"gapFreeConsensus\": \"" << cs << "\"," << std::endl;
      rfile << "\"gappedConsensus\": \"" << gapped << "\"," << std::endl;
      rfile << "\"msa\": " << std::endl;
      rfile << "[" << std::endl;
      for(uint32_t i = 0; i < (TAIndex) align.shape()[0]; ++i) {
	if (i!=0) rfile << ',' << std::endl;
	alignedTraceByRow(rfile, align, i, c.ab[idxMap[seqidx[i]]].stem().string(), fwd[seqidx[i]], false);
      }
      rfile << "]," << std::endl;
      rfile << "\"gappedTraces\": " << std::endl;
      rfile << "[" << std::endl;
      for(uint32_t i = 0; i < (TAIndex) align.shape()[0]; ++i) {
	if (i!=0) rfile << ", ";
	Trace tr;
	int32_t ft = traceFormat(c.ab[idxMap[seqidx[i]]].string());
	if (ft == 0) {
	  if (!readab(c.ab[idxMap[seqidx[i]]].string(), tr)) return -1;
	} else if (ft == 1) {
	    if (!readscf(c.ab[idxMap[seqidx[i]]].string(), tr)) return -1;
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
	if (c.trimStringency) trimTrace(c, bc, trimLeft, trimRight);
	
	// Hard trim of the trace data structures
	BaseCalls nbc;
	trimTrace(tr, bc, trimLeft, trimRight, nbc);

	// Reverse complement trace if necessary
	BaseCalls padbc;
	Trace padtr;
	if (fwd[seqidx[i]]) alignmentTracePadding(align, tr, nbc, i, padtr, padbc);
	else {
	  BaseCalls tbc;
	  Trace ttr;
	  // Reverese complement
	  reverseComplementTrace(tr, nbc, ttr, tbc);
	  alignmentTracePadding(align, ttr, tbc, i, padtr, padbc);
	}
	// Append gapped trace to output
	assemblyTrace(rfile, padbc, padtr, c.ab[idxMap[seqidx[i]]].stem().string());
      }
      rfile << "]" << std::endl;
      rfile << "}" << std::endl;
      rfile.close();
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

    // Output consensus
    if (c.format == "fasta") {
      filename = c.outprefix + ".cons.fa";
      std::ofstream csfile(filename.c_str());
      csfile << ">Consensus" << std::endl;
      csfile << cs << std::endl;
      csfile.close();
    } else if (c.format == "fastq") {
      filename = c.outprefix + ".cons.fq";
      std::ofstream csfile(filename.c_str());
      csfile << "@Consensus" << std::endl;
      csfile << cs << std::endl;
      csfile << "+" << std::endl;
      csfile << qstr << std::endl;
      csfile.close();
    }

    // Done
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Done." << std::endl;
    return 0;
  }

}

#endif
