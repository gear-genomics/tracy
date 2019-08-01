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
#include "fasta.h"
#include "trim.h"
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
    int32_t gapopen;
    int32_t gapext;
    int32_t match;
    int32_t mismatch;
    float pratio;
    float trimStringency;
    std::string outprefix;
    DnaScore<int32_t> aliscore;
    boost::filesystem::path align;
    boost::filesystem::path ab;
    boost::filesystem::path outfile;
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
      ("reference,r", boost::program_options::value<boost::filesystem::path>(&c.genome), "(gzipped) fasta or wildtype ab1 file")
      ("pratio,p", boost::program_options::value<float>(&c.pratio)->default_value(0.33), "peak ratio to call base")
      ("kmer,k", boost::program_options::value<uint16_t>(&c.kmer)->default_value(15), "kmer size to anchor trace")
      ("support,s",  boost::program_options::value<uint16_t>(&c.minKmerSupport)->default_value(3), "min. kmer support")
      ("maxindel,i", boost::program_options::value<uint16_t>(&c.maxindel)->default_value(1000), "max. indel size in Sanger trace")
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
      ("trimLeft,q", boost::program_options::value<uint16_t>(&c.trimLeft)->default_value(50), "trim size left")
      ("trimRight,u", boost::program_options::value<uint16_t>(&c.trimRight)->default_value(50), "trim size right")
      ;

    boost::program_options::options_description otp("Output options");
    otp.add_options()
      ("linelimit,l", boost::program_options::value<uint16_t>(&c.linelimit)->default_value(60), "alignment line length")
      ("outprefix,o", boost::program_options::value<std::string>(&c.outprefix)->default_value("out"), "output prefix")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.ab), "ab1")
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
      std::cout << "Usage: tracy " << argv[0] << " [OPTIONS] -r genome.fa trace.ab1" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }
    if (c.maxindel < 1) c.maxindel = 1;

    // Check trimming parameters
    if (c.trimStringency > 9) c.trimStringency = 9;
    
    // Check reference
    if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
      std::cerr << "Reference file is missing: " << c.genome.filename().string() << std::endl;
      return 1;
    }

    // Check input trace
    if (!(boost::filesystem::exists(c.ab) && boost::filesystem::is_regular_file(c.ab) && boost::filesystem::file_size(c.ab))) {
      std::cerr << "Input trace file is missing: " << c.ab.filename().string() << std::endl;
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

    // Any basecalls
    if (!tr.basecallpos.size()) {
      std::cerr << "Trace file lacks basecalls!" << std::endl;
      return -1;
    }

    // Alignment options
    AlignConfig<true, false> semiglobal;
    c.aliscore = DnaScore<int32_t>(c.match, c.mismatch, c.gapopen, c.gapext);

    // Call bases
    BaseCalls bc;
    basecall(tr, bc, c.pratio);

    // Get trim sizes
    if (c.trimStringency >= 1) {
      uint32_t trimLeft = 0;
      uint32_t trimRight = 0;
      trimTrace(c, bc, trimLeft, trimRight);
      c.trimLeft = trimLeft;
      c.trimRight = trimRight;
    }

    // Check trim sizes
    if (c.trimLeft + c.trimRight >= bc.bcPos.size()) {
      std::cerr << "The sum of the left and right trim size is larger than the trace!" << std::endl;
      return -1;
    }
    
    // Output trace information
    traceTxtOut(c.outprefix + ".abif", bc, tr, c.trimLeft, c.trimRight);
    
    // Full trace profile
    typedef boost::multi_array<float, 2> TProfile;
    TProfile fulltraceprofile;
    createProfile(tr, bc, fulltraceprofile);

    // Create trimmed trace profile
    TProfile trimmedtrace;
    createProfile(tr, bc, trimmedtrace, c.trimLeft, c.trimRight);

    // Load reference
    ReferenceSlice rs;
    TProfile referenceprofile;
    rs.filetype = genomeType(c.genome.string());
    if (rs.filetype == -1) {
      std::cerr << "Unknown reference file format!" << std::endl;
      return -1;
    }

    // Find reference match
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Find reference match" << std::endl;
    if ((rs.filetype == 0) || (rs.filetype == 1)) {
      // Preliminary alignment
      typedef boost::multi_array<char, 2> TAlign;
      TAlign align;
      TProfile prefslice;	  
      if (rs.filetype == 0) {
	// Indexed genome
	csa_wt<> fm_index;
	if (!loadFMIdx(c, rs, fm_index)) return -1;
	if (!getReferenceSlice(c, fm_index, bc, rs)) return -1;
	createProfile(c, rs, prefslice);
      } else {
	// Single FASTA
	std::string faname = "";
	std::string seq = "";
	if (!loadSingleFasta(c.genome.string(), faname, seq)) return -1;
	if (seq.size() > MAX_SINGLE_FASTA_SIZE) {
	  std::cerr << "Reference is larger than 50Kbp. Please use a smaller reference slice or an indexed genome!" << std::endl;
	  return -1;
	}

	// Profile
	TProfile fwdprofile;
	_createProfile(seq, fwdprofile);
	TProfile revprofile;
	reverseComplementProfile(fwdprofile, revprofile);
	
	// Alignment scores
	int32_t gsFwd = gotohScore(trimmedtrace, fwdprofile, semiglobal, c.aliscore);
	int32_t gsRev = gotohScore(trimmedtrace, revprofile, semiglobal, c.aliscore);
	
	// Forward or reverse?
	rs.kmersupport = 0;
	rs.pos = 0;
	rs.chr = faname;
	rs.refslice = seq;
	if (gsFwd > gsRev) {
	  rs.forward = true;
	  copyProfile(fwdprofile, prefslice);
	} else {
	  rs.forward = false;
	  reverseComplement(rs.refslice);
	  copyProfile(revprofile, prefslice);
	}
      }
      
      // Align trimmed trace to profile
      gotoh(trimmedtrace, prefslice, align, semiglobal, c.aliscore);
      trimReferenceSlice(c, align, rs);
      createProfile(c, rs, referenceprofile);
    } else if (rs.filetype == 2) {
      // Wildtype trace
      Trace gtr;
      int32_t gft = traceFormat(c.genome.string());
      if (gft == 0) {
	if (!readab(c.genome.string(), gtr)) return -1;
      } else if (gft == 1) {
	if (!readscf(c.genome.string(), gtr)) return -1;
      } else {
	std::cerr << "Unknown trace file type!" << std::endl;
	return -1;
      }

      // Basecalling and reference profile
      BaseCalls gbc;
      basecall(gtr, gbc, c.pratio);

      // Figure out if fwd or rev
      TProfile fwdprofile;
      createProfile(gtr, gbc, fwdprofile);
      TProfile revprofile;
      reverseComplementProfile(fwdprofile, revprofile);

      // Alignment scores
      int32_t gsFwd = gotohScore(trimmedtrace, fwdprofile, semiglobal, c.aliscore);
      int32_t gsRev = gotohScore(trimmedtrace, revprofile, semiglobal, c.aliscore);
      
      // Forward or reverse?
      rs.kmersupport = 0;
      rs.pos = 0;
      rs.chr = "wildtype";
      rs.refslice = gbc.primary;
      if (gsFwd > gsRev) {
	rs.forward = true;
	copyProfile(fwdprofile, referenceprofile);
      } else {
	rs.forward = false;
	reverseComplement(rs.refslice);
	copyProfile(revprofile, referenceprofile);
      }
    } else {
      std::cerr << "Unknown reference file type!" << std::endl;
      return - 1;
    }

    // Semi-global alignment
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Alignment" << std::endl;
    typedef boost::multi_array<char, 2> TAlign;
    TAlign final;
    int32_t score = gotoh(fulltraceprofile, referenceprofile, final, semiglobal, c.aliscore);
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

    // Horizontal alignment
    std::string alignfilename = c.outprefix + ".align.fa";
    std::ofstream vfile(alignfilename.c_str());
    typedef typename TAlign::index TAIndex;
    vfile << ">" << c.ab.stem().string() << std::endl;
    for(TAIndex j = 0; j < (TAIndex) final.shape()[1]; ++j) vfile << final[0][j];
    vfile << std::endl;
    vfile << ">" << rs.chr;
    if (rs.forward) vfile << " (forward)" << std::endl;
    else vfile << " (reverse)" << std::endl;
    for(TAIndex j = 0; j < (TAIndex) final.shape()[1]; ++j) vfile << final[1][j];
    vfile << std::endl;
    vfile.close();

    // Show ClustalW like alignment
    plotAlignment(c.outprefix + ".txt", final, rs, score, c.linelimit);

    // JSON output
    traceAlignJsonOut(c.outprefix + ".json", nbc, ntr, rs, final);

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


