#ifndef CONSENSUS_H
#define CONSENSUS_H

#define BOOST_DISABLE_ASSERTS
#include <boost/math/special_functions/round.hpp>
#include <boost/math/distributions/normal.hpp>
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

#include <sdsl/suffix_arrays.hpp>

#include <htslib/faidx.h>

#include "abif.h"
#include "align.h"
#include "gotoh.h"
#include "fasta.h"
#include "trim.h"
#include "fmindex.h"
#include "json.h"
#include "profile.h"

using namespace sdsl;

namespace tracy {

#define SMALLEST_GL -1000

  struct ConsensusConfig {
    bool computeUnion;
    bool useIUPAC;
    uint16_t linelimit;
    uint16_t trimLeft1;
    uint16_t trimRight1;
    uint16_t trimLeft2;
    uint16_t trimRight2;
    int32_t gapopen;
    int32_t gapext;
    int32_t match;
    int32_t mismatch;
    float pratio;
    float trimStringency;
    std::string label;
    std::string outprefix;
    DnaScore<int32_t> aliscore;
    boost::filesystem::path align;
    boost::filesystem::path outfile;
    std::vector<boost::filesystem::path> files;
  };

  template<typename TConfig>
  inline void
  consensusFastaOut(TConfig const& c, std::string const& cons) {
    // Output trace
    std::string outfile = c.outprefix + ".fa";
    std::ofstream rfile(outfile.c_str());
    rfile << ">" << c.label << std::endl;
    rfile << cons << std::endl;
    rfile.close();  
  }

  template<typename TConfig>
  inline void
  consensusFastqOut(TConfig const& c, std::string const& cons, std::vector<uint32_t>& qual) {
    // Output trace
    std::string outfile = c.outprefix + ".fq";
    std::ofstream rfile(outfile.c_str());
    rfile << "@" << c.label << std::endl;
    rfile << cons << std::endl;
    rfile << "+" << std::endl;
    for(uint32_t i = 0; i < qual.size(); ++i) {
      int32_t qval = qual[i] + 33;
      if (qval > 122) qval = 122;
      rfile << (char) (qval);
    }
    rfile << std::endl;
    rfile.close();
  }

  template<typename TConfig>
  inline void
  gtLetter(TConfig const& c, std::vector<double>& cl, std::string& cons, std::vector<uint32_t>& qual) {
    // Genotype likelihoods
    std::vector<double> gl(6);
    double total = 0;
    for(uint32_t k = 0; k < cl.size(); ++k) total += cl[k];
    for(uint32_t k = 0; k < cl.size(); ++k) {
      if (total > 0) cl[k] /= total;
      else cl[k] = 0;
      if (cl[k] > 0) {
	gl[k] = std::log10(cl[k]);
	if (gl[k] < SMALLEST_GL) gl[k] = SMALLEST_GL;
      } else gl[k] = SMALLEST_GL;
    }

    // Rescale by best consensus letter
    uint32_t glBest = 0;
    uint32_t gl2ndBest = 1;
    if (gl[glBest] < gl[gl2ndBest]) {
      glBest = 1;
      gl2ndBest = 0;
    }
    for(uint32_t k = 2; k < gl.size(); ++k) {
      if (gl[k] > gl[glBest]) {
	gl2ndBest = glBest;
	glBest = k;
      } else if (gl[k] > gl[gl2ndBest]) {
	gl2ndBest = k;
      }
    }
    double glBestVal = gl[glBest];
    bool ambiguous = false;
    if (c.useIUPAC) {
      if (gl[gl2ndBest] > -1) {
	// Make sure both nucleotides are A, C, G, T
	if ((glBest <= 3) && (gl2ndBest <= 3)) ambiguous = true;
      }
    }
    for(uint32_t k = 0; k < gl.size(); ++k) gl[k] -= glBestVal;

    // Compute quality (2nd best basecall to best basecall)
    uint32_t bestPL = (uint32_t) boost::math::round(-10 * gl[glBest]);
    uint32_t best2ndPL = (uint32_t) boost::math::round(-10 * gl[gl2ndBest]);
    double likelihood = std::log10(1 - 1 / (std::pow((double) 10, -( (double) bestPL / (double) 10)) + std::pow((double) 10, -( (double) best2ndPL / (double) 10))));
    likelihood = (likelihood > SMALLEST_GL) ? likelihood : SMALLEST_GL;
    int32_t gqval = (int32_t) boost::math::round(-10 * likelihood);
    if (gqval < 0) gqval = 0;

    // Debug
    //for(uint32_t k = 0; k < gl.size(); ++k) std::cerr << k << ',' << gl[k] << std::endl;
    //std::cerr << gqval << std::endl;

    // Determine consensus letter and quality
    // 'A', 'C', 'G', 'T', 'N', '-'
    if (ambiguous) {
      // Use IUPAC
      char c1 = 'A';
      if (glBest == 0) c1 = 'A';
      else if (glBest == 1) c1 = 'C';
      else if (glBest == 2) c1 = 'G';
      else if (glBest == 3) c1 = 'T';
      char c2 = 'A';
      if (gl2ndBest == 0) c2 = 'A';
      else if (gl2ndBest == 1) c2 = 'C';
      else if (gl2ndBest == 2) c2 = 'G';
      else if (gl2ndBest == 3) c2 = 'T';
      cons += iupac(c1, c2);
    } else {
      if (glBest == 0) cons += 'A';
      else if (glBest == 1) cons += 'C';
      else if (glBest == 2) cons += 'G';
      else if (glBest == 3) cons += 'T';
      else if (glBest == 4) cons += 'N';
      else cons += '-';
    }
    qual.push_back(gqval);
  }
  
  template<typename TConfig, typename TProfile>
  inline void
  consLetter(TConfig const& c, TProfile const& p1, TProfile const& p2, int32_t const s1, int32_t const s2, std::string& cons, std::vector<uint32_t>& qual) {
    std::vector<double> cl(6);
    for(uint32_t k = 0; k < 6; ++k) cl[k] = p1[k][s1] + p2[k][s2];
    gtLetter(c, cl, cons, qual); 
  }

  template<typename TConfig, typename TProfile>
  inline void
  consLetter(TConfig const& c, TProfile const& p, int32_t const s, std::string& cons, std::vector<uint32_t>& qual) {
    std::vector<double> cl(6);
    for(uint32_t k = 0; k < 6; ++k) cl[k] = p[k][s];
    gtLetter(c, cl, cons, qual); 
  }
  
  template<typename TConfig, typename TAlign, typename TProfile>
  inline void
  pairwiseConsensus(TConfig const& c, TAlign const& align, TProfile const& trimmedtrace1, TProfile const& trimmedtrace2, std::string& cons, std::vector<uint32_t>& qual) {
    typedef typename TAlign::index TAIndex;
    // Parse alignment
    int32_t seq1 = 0;
    int32_t seq2 = 0;
    uint32_t gapEx = 0;
    uint32_t mm = 0; // Number of mismatches
    uint32_t ma = 0; // Number of matches
    uint32_t go = 0; // Number of gap-openings
    uint32_t ge = 0; // Number of gap-extensions
    bool inGap=false;
    for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
      if ((align[0][j] == '-') || (align[1][j] == '-')) {
	if ((seq1) && (seq2)) {
	  if (!inGap) {
	    inGap = true;
	    gapEx = 0;
	  }
	  gapEx += 1;
	}
	if (align[0][j] != '-') {
	  if (c.computeUnion) consLetter(c, trimmedtrace1, seq1, cons, qual);
	  ++seq1;
	}
	if (align[1][j] != '-') {
	  if (c.computeUnion) consLetter(c, trimmedtrace2, seq2, cons, qual);
	  ++seq2;
	}
      } else {
	if (inGap) {
	  ge += gapEx;
	  go += 1;
	  inGap=false;
	}
	if (align[0][j] == align[1][j]) ma += 1;
	else mm += 1;
	//std::cerr << align[0][j] << ',' << align[1][j] << std::endl;
	//for(uint32_t k = 0; k < 6; ++k) std::cerr << k << ':' << trimmedtrace1[k][seq1] << ',' << trimmedtrace2[k][seq2] << std::endl;
	consLetter(c, trimmedtrace1, trimmedtrace2, seq1, seq2, cons, qual);
	++seq1;
	++seq2;
      }
    }
    // Debug
    //std::cerr << go << ',' << ge << ',' << ma << ',' << mm << std::endl;
    //std::cerr << seq1 << ',' << trimmedtrace1.shape()[1] << std::endl;
    //std::cerr << seq2 << ',' << trimmedtrace2.shape()[1] << std::endl;
  }
  

  template<typename TConfig, typename TAlign>
  inline void
  plotClustalPairwise(TConfig const& c, TAlign const& align, bool const forward, int32_t const score, uint32_t const linelimit) {
    typedef typename TAlign::index TAIndex;

    // Sequence header
    uint32_t fald = linelimit + 14;
    std::string filename = c.outprefix + ".txt";
    std::ofstream ofile(c.outprefix + ".txt");
    ofile << ">" << c.files[0].stem().string() << std::endl;
    int32_t count = 0;
    for(TAIndex j = 0; j< (TAIndex) align.shape()[1]; ++j) {
      if (align[0][j] != '-')  {
	ofile << align[0][j];
	if ((count+1) % fald == 0) ofile << std::endl;
	++count;
      }
    }
    if (count % fald != 0) ofile << std::endl;
    ofile << ">" << c.files[1].stem().string();
    if (forward) ofile << " (forward)" << std::endl;
    else ofile << " (reverse)" << std::endl;
    count = 0;
    for(TAIndex j = 0; j< (TAIndex) align.shape()[1]; ++j) {
      if (align[1][j] != '-') {
	ofile << align[1][j];
	if ((count+1)% fald == 0) ofile << std::endl;
	++count;
      }
    }
    if (count % fald != 0) ofile << std::endl;
    ofile << std::endl;
    ofile << "Alignment score: " << score << std::endl;
    ofile << "#";
    for(uint32_t i = 1; i < fald; ++i) ofile << "-";
    ofile << std::endl;
    ofile << std::endl;

    // Alignment
    int32_t vi = 1;
    int32_t ri = 1;
    std::string f1 = c.files[0].stem().string().substr(0,8);
    while (f1.size() < 8) f1 += " ";
    std::string f2 = c.files[1].stem().string().substr(0,8);
    while (f2.size() < 8) f2 += " ";
    uint32_t blockcount = 0;
    int32_t s = 0;
    int32_t e = (TAIndex) align.shape()[1];
    while (s < e) {
      ofile << f1 << std::setw(5) << vi << ' ';
      for(TAIndex j = s; ((j < (TAIndex) e) && (j < s + linelimit)); ++j) {
	ofile << align[0][j];
	if (align[0][j] != '-') ++vi;
      }
      ofile << std::endl;
      ofile << "              ";
      for(TAIndex j = s; ((j < (TAIndex) e) && (j < s + linelimit)); ++j) {
	if (align[0][j] == align[1][j]) ofile << "|";
	else ofile << " ";
      }
      ofile << std::endl;
      ofile << f2 << std::setw(5) << ri << ' ';
      for(TAIndex j = s; ((j < (TAIndex) e) && (j < s + linelimit)); ++j) {
	ofile << align[1][j];
	if (align[1][j] != '-') ++ri;
      }
      ofile << std::endl;
      ofile << std::endl;
      s += linelimit;
      ++blockcount;
    }
    if (blockcount < 6) {
      // Add spacer for small alignments
      for(uint32_t i = blockcount; i < 6; ++i) {
	for(uint32_t k = 0; k<4; ++k) ofile << std::endl;
      }
    }
    ofile << "#";
    for(uint32_t i = 1; i < fald; ++i) ofile << "-";
    ofile << std::endl;
    ofile << "#";
    for(uint32_t i = 1; i < fald; ++i) ofile << "-";
    ofile << std::endl;
    ofile << std::endl;
    ofile << std::endl;
    ofile.close();
  }



  
  int consensus(int argc, char** argv) {

#ifdef PROFILE
    ProfilerStart("tracy.prof");
#endif
    
    ConsensusConfig c;
  
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("label,b", boost::program_options::value<std::string>(&c.label)->default_value("Consensus"), "sample label")
      ("pratio,p", boost::program_options::value<float>(&c.pratio)->default_value(0.33), "peak ratio to call base")
      ;

    boost::program_options::options_description alignment("Alignment options");
    alignment.add_options()
      ("gapopen,g", boost::program_options::value<int32_t>(&c.gapopen)->default_value(-10), "gap open")
      ("gapext,e", boost::program_options::value<int32_t>(&c.gapext)->default_value(-4), "gap extension")
      ("match,m", boost::program_options::value<int32_t>(&c.match)->default_value(3), "match")
      ("mismatch,n", boost::program_options::value<int32_t>(&c.mismatch)->default_value(-5), "mismatch")
      ;

    boost::program_options::options_description tro("Trimming options");
    tro.add_options()
      ("trim,t", boost::program_options::value<float>(&c.trimStringency)->default_value(0), "trimming stringency [1:9], 0: use trimLeft and trimRight")
      ("trimLeft1,q", boost::program_options::value<uint16_t>(&c.trimLeft1)->default_value(50), "trim size left (1st trace)")
      ("trimRight1,u", boost::program_options::value<uint16_t>(&c.trimRight1)->default_value(50), "trim size right (1st trace)")
      ("trimLeft2,r", boost::program_options::value<uint16_t>(&c.trimLeft2)->default_value(50), "trim size left (2nd trace)")
      ("trimRight2,s", boost::program_options::value<uint16_t>(&c.trimRight2)->default_value(50), "trim size right (2nd trace)")
      ;

    boost::program_options::options_description otp("Output options");
    otp.add_options()
      ("linelimit,l", boost::program_options::value<uint16_t>(&c.linelimit)->default_value(60), "alignment line length")
      ("outprefix,o", boost::program_options::value<std::string>(&c.outprefix)->default_value("out"), "output prefix")
      ("intersect,i", "use only trace intersection for consensus")
      ("iupac,a", "use IUPAC nucleotide code in consensus (max. 2 nucleotides)")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<std::vector<boost::filesystem::path> >(&c.files), "ab1")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(alignment).add(tro).add(otp).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(alignment).add(tro).add(otp);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    
    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file"))) {
      std::cout << "Usage: tracy " << argv[0] << " [OPTIONS] trace1.ab1 trace2.ab1" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }

    // 2 files?
    if (c.files.size() != 2) {
      std::cerr << "Exactly 2 input trace files are required!" << std::endl;
      return 1;
    }
    
    // Check input trace
    for(uint32_t i = 0; i < c.files.size(); ++i) {
      if (!(boost::filesystem::exists(c.files[i]) && boost::filesystem::is_regular_file(c.files[i]) && boost::filesystem::file_size(c.files[i]))) {
	std::cerr << "Input trace file is missing: " << c.files[i].filename().string() << std::endl;
	return 1;
      }
    }

    // Intersection for consensus
    if (vm.count("intersect")) c.computeUnion = false;
    else c.computeUnion = true;

    // Use IUPAC
    if (vm.count("iupac")) c.useIUPAC = true;
    else c.useIUPAC = false;
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "tracy ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    // Load first trace file
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Load " << c.files[0].string() << " file" << std::endl;
    Trace tr1;
    int32_t ft1 = traceFormat(c.files[0].string());
    if (ft1 == 0) {
      if (!readab(c.files[0].string(), tr1)) return -1;
    } else if (ft1 == 1) {
      if (!readscf(c.files[0].string(), tr1)) return -1;
    } else {
      std::cerr << "Unknown trace file type!" << std::endl;
      return -1;
    }
    // Any basecalls
    if (!tr1.basecallpos.size()) {
      std::cerr << "Trace file lacks basecalls!" << std::endl;
      return -1;
    }
    
    // Load second trace file
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Load " << c.files[1].string() << " file" << std::endl;
    Trace tr2;
    int32_t ft2 = traceFormat(c.files[1].string());
    if (ft2 == 0) {
      if (!readab(c.files[1].string(), tr2)) return -1;
    } else if (ft2 == 1) {
      if (!readscf(c.files[1].string(), tr2)) return -1;
    } else {
      std::cerr << "Unknown trace file type!" << std::endl;
      return -1;
    }
    // Any basecalls
    if (!tr2.basecallpos.size()) {
      std::cerr << "Trace file lacks basecalls!" << std::endl;
      return -1;
    }

    // Alignment options
    AlignConfig<true, true> global;
    c.aliscore = DnaScore<int32_t>(c.match, c.mismatch, c.gapopen, c.gapext);

    // Call bases
    BaseCalls bc1;
    basecall(tr1, bc1, c.pratio);
    BaseCalls bc2;
    basecall(tr2, bc2, c.pratio);

    // Trace trimming
    if (c.trimStringency >= 1) {
      uint32_t trimLeft = 0;
      uint32_t trimRight = 0;
      trimTrace(c, bc1, trimLeft, trimRight);
      c.trimLeft1 = trimLeft;
      c.trimRight1 = trimRight;
      trimLeft = 0;
      trimRight = 0;
      trimTrace(c, bc2, trimLeft, trimRight);
      c.trimLeft2 = trimLeft;
      c.trimRight2 = trimRight;
    }

    // Check trim sizes
    if (c.trimLeft1 + c.trimRight1 >= bc1.bcPos.size()) {
      std::cerr << "The sum of the left and right trim size is larger than the trace: " << c.files[0].string() << std::endl;
      return -1;
    }
    if (c.trimLeft2 + c.trimRight2 >= bc2.bcPos.size()) {
      std::cerr << "The sum of the left and right trim size is larger than the trace:" << c.files[1].string() << std::endl;
      return -1;
    }
    
    // Output trace information
    traceTxtOut(c.outprefix + "_1st.abif", bc1, tr1, c.trimLeft1, c.trimRight1);
    traceTxtOut(c.outprefix + "_2nd.abif", bc2, tr2, c.trimLeft2, c.trimRight2);
    
    // Full trace profile
    typedef boost::multi_array<float, 2> TProfile;
    TProfile fulltraceprofile1;
    createProfile(tr1, bc1, fulltraceprofile1);
    TProfile fulltraceprofile2;
    createProfile(tr2, bc2, fulltraceprofile2);

    // Create trimmed trace profile
    TProfile trimmedtrace1;
    createProfile(tr1, bc1, trimmedtrace1, c.trimLeft1, c.trimRight1);
    TProfile fwdprofile2;
    createProfile(tr2, bc2, fwdprofile2, c.trimLeft2, c.trimRight2);
    TProfile revprofile2;
    reverseComplementProfile(fwdprofile2, revprofile2);

    // Alignment scores
    int32_t gsFwd = gotohScore(trimmedtrace1, fwdprofile2, global, c.aliscore);
    int32_t gsRev = gotohScore(trimmedtrace1, revprofile2, global, c.aliscore);
      
    // Forward or reverse?
    bool forward = true;
    TProfile trimmedtrace2;
    if (gsFwd > gsRev) {
      copyProfile(fwdprofile2, trimmedtrace2);
    } else {
      forward = false;
      copyProfile(revprofile2, trimmedtrace2);
    }

    // Semi-global alignment
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Alignment" << std::endl;
    typedef boost::multi_array<char, 2> TAlign;
    TAlign fali;
    int32_t score = gotoh(trimmedtrace1, trimmedtrace2, fali, global, c.aliscore);
    // Debug Alignment
    //for(uint32_t i = 0; i < fali.shape()[0]; ++i) {
    //for(uint32_t j = 0; j < fali.shape()[1]; ++j) std::cerr << fali[i][j];
    //std::cerr << std::endl;
    //}    

    // Output
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] " << "Output" << std::endl;

    // Horizontal alignment
    std::string alignfilename = c.outprefix + ".align.fa";
    std::ofstream vfile(alignfilename.c_str());
    typedef typename TAlign::index TAIndex;
    vfile << ">" << c.files[0].stem().string() << std::endl;
    for(TAIndex j = 0; j < (TAIndex) fali.shape()[1]; ++j) vfile << fali[0][j];
    vfile << std::endl;
    vfile << ">" << c.files[1].stem().string();
    if (forward) vfile << " (forward)" << std::endl;
    else vfile << " (reverse)" << std::endl;
    for(TAIndex j = 0; j < (TAIndex) fali.shape()[1]; ++j) vfile << fali[1][j];
    vfile << std::endl;
    vfile.close();

    // Consensus
    std::string cons;
    std::vector<uint32_t> qual;
    pairwiseConsensus(c, fali, trimmedtrace1, trimmedtrace2, cons, qual);

    // Output
    consensusFastaOut(c, cons);
    consensusFastqOut(c, cons, qual);

    // Show ClustalW like alignment
    plotClustalPairwise(c, fali, forward, score, c.linelimit);

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


