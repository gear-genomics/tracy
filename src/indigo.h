#ifndef INDIGO_H
#define INDIGO_H

#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#include "decompose.h"
#include "trim.h"
#include "web.h"
#include "variants.h"
#include "fmindex.h"

using namespace sdsl;

namespace tracy {
  
  struct IndigoConfig {
    bool callvariants;
    bool annotatevariants;
    uint16_t linelimit;
    uint16_t trimLeft;
    uint16_t trimRight;
    uint16_t kmer;
    uint16_t maxindel;
    uint16_t madc;
    uint16_t minKmerSupport;
    uint16_t qualCut;
    int32_t gapopen;
    int32_t gapext;
    int32_t match;
    int32_t mismatch;
    float pratio;
    float trimStringency;
    std::string annotate;
    std::string host;
    std::string outprefix;
    DnaScore<int32_t> aliscore;
    boost::filesystem::path outfile;
    boost::filesystem::path ab;
    boost::filesystem::path genome;
  };

  int indigo(int argc, char** argv) {
    IndigoConfig c;
  
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("genome,r", boost::program_options::value<boost::filesystem::path>(&c.genome), "(gzipped) fasta or wildtype ab1 file")
      ("pratio,p", boost::program_options::value<float>(&c.pratio)->default_value(0.33), "peak ratio to call base")
      ("kmer,k", boost::program_options::value<uint16_t>(&c.kmer)->default_value(15), "kmer size")
      ("support,s",  boost::program_options::value<uint16_t>(&c.minKmerSupport)->default_value(3), "min. kmer support")
      ("maxindel,i", boost::program_options::value<uint16_t>(&c.maxindel)->default_value(1000), "max. indel size in Sanger trace")
      ("annotate,a", boost::program_options::value<std::string>(&c.annotate), "annotate variants [homo_sapiens|homo_sapiens_hg19|mus_musculus|danio_rerio|...]")
      ("callVariants,v", "call variants in trace")
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
      ("madc,c", boost::program_options::value<uint16_t>(&c.madc)->default_value(5), "MAD cutoff")
      ("qualCut,z",  boost::program_options::value<uint16_t>(&c.qualCut)->default_value(45), "variant calling quality threshold")
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
      std::cout << "Usage: tracy " << argv[0] << " [OPTIONS] trace.ab1" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }
    if (c.maxindel < 1) c.maxindel = 1;

    // Check trimming parameters
    if (c.trimStringency > 9) c.trimStringency = 9;
    
    // Variant calling
    if (vm.count("callVariants")) c.callvariants = true;
    else c.callvariants = false;
    if (vm.count("annotate")) {
      c.annotatevariants = true;
      c.callvariants = true;
      c.host = "rest.ensembl.org";
      c.annotate = fixSpeciesName(c.annotate);
      // hg19 workaround
      if (c.annotate == "homo_sapiens_hg19") {
	c.host = "grch37.rest.ensembl.org";
	c.annotate = "homo_sapiens";
      }
      if (!speciesExist(c.annotate)) c.annotatevariants = false;
    } else c.annotatevariants = false;
    
    // Check ab1
    if (!(boost::filesystem::exists(c.ab) && boost::filesystem::is_regular_file(c.ab) && boost::filesystem::file_size(c.ab))) {
      std::cerr << "Trace file is missing: " << c.ab.string() << std::endl;
      return 1;
    }
    
    // Check reference
    if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
      std::cerr << "Reference file is missing: " << c.genome.filename().string() << std::endl;
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

    // Create trimmed trace profile
    typedef boost::multi_array<float, 2> TProfile;
    TProfile trimmedtrace;
    createProfile(tr, bc, trimmedtrace, c.trimLeft, c.trimRight);

    // Identify position of indel shift in Sanger trace
    TraceBreakpoint bp;
    findBreakpoint(trimmedtrace, bp);

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
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Find Reference Match" << std::endl;
    TProfile prefslice;
    if ((rs.filetype == 0) || (rs.filetype == 1)) {
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
	copyProfile(fwdprofile, prefslice);
      } else {
	rs.forward = false;
	reverseComplement(rs.refslice);
	copyProfile(revprofile, prefslice);
      }
    } else {
      std::cerr << "Unknown reference file type!" << std::endl;
      return - 1;
    }

    // Align trimmed trace to profile    
    typedef boost::multi_array<char, 2> TAlign;
    TAlign align;
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Alignment" << std::endl;
    int32_t aliTrimScore = gotoh(trimmedtrace, prefslice, align, semiglobal, c.aliscore);
    double seqsize = trimmedtrace.shape()[1];
    double matchFraction = 0.35;
    double scoreThreshold = seqsize * matchFraction * c.aliscore.match + seqsize * (1 - matchFraction) * c.aliscore.mismatch;
    if (aliTrimScore <= scoreThreshold) {
      std::cerr << "Alignment of trace to reference failed!" << std::endl;
      return -1;
    }
    
    // Hom. InDel search
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "InDel Search" << std::endl;
    if (!bp.indelshift) {
      // Find breakpoint for hom. indels
      if (!findHomozygousBreakpoint(align, bp)) return -1;
    }

    // Debug Breakpoint & Alignment
    //std::cerr << "Breakpoint: " << bp.indelshift << ',' << bp.traceleft << ',' << bp.breakpoint << ',' << bp.bestDiff << std::endl;
    //for(uint32_t i = 0; i<align.shape()[0]; ++i) {
    //uint32_t alignedNuc = 0;
    //for(uint32_t j = 0; j<align.shape()[1]; ++j) {
    //if (align[0][j] != '-') {
    //++alignedNuc;
    //if (alignedNuc == bp.breakpoint) std::cerr << "#####";
    //}
    //std::cerr << align[i][j];
    //}
    //std::cerr << std::endl;
    //}    

    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Decompose Chromatogram" << std::endl;
    
    // Decompose alleles
    typedef std::pair<int32_t, int32_t> TIndelError;
    typedef std::vector<TIndelError> TDecomposition;
    TDecomposition dcp;
    if (!decomposeAlleles(c, align, bc, bp, rs, dcp)) return -1;
    writeDecomposition(c.outprefix + ".decomp", dcp);

    // Generate plain nucleotide sequence for second allele
    generateSecondaryDecomposed(tr, bc);

    // Estimate allelic fractions
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Estimate allelic fractions" << std::endl;
    typedef std::pair<double, double> TFractions;
    TFractions a1a2 = allelicFraction(c, tr, bc);
    
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Allele-specific alignments" << std::endl;
    
    // Allele1
    typedef boost::multi_array<char, 2> TAlign;
    TAlign alignPrimary;
    std::string pri = trimmedSeq(bc.primary, c.trimLeft, c.trimRight);
    gotoh(pri, rs.refslice, alignPrimary, semiglobal, c.aliscore);
    // Trim initial reference slice
    ReferenceSlice allele1(rs);
    trimReferenceSlice(c, alignPrimary, allele1);
    typedef boost::multi_array<char, 2> TAlign;
    TAlign final1;
    int32_t a1Score = gotoh(pri, allele1.refslice, final1, semiglobal, c.aliscore);
    plotAlignment(c.outprefix + ".align1", final1, allele1, 1, a1Score, a1a2, c.linelimit);

    // Allele2
    TAlign alignSecondary;
    std::string sec = trimmedSeq(bc.secDecompose, c.trimLeft, c.trimRight);
    gotoh(sec, rs.refslice, alignSecondary, semiglobal, c.aliscore);
    // Trim initial reference slice
    ReferenceSlice allele2(rs);
    trimReferenceSlice(c, alignSecondary, allele2);
    TAlign final2;
    int32_t a2Score = gotoh(sec, allele2.refslice, final2, semiglobal, c.aliscore);
    plotAlignment(c.outprefix + ".align2", final2, allele2, 2, a2Score, a1a2, c.linelimit);

    // Allele1 vs. Allele2
    TAlign final3;
    AlignConfig<false, false> global;
    ReferenceSlice secrs;
    secrs.refslice = sec;
    secrs.forward = 1;
    secrs.pos = 0;
    secrs.chr = "Alt2";
    int32_t a3Score = gotoh(pri, secrs.refslice, final3, global, c.aliscore);
    plotAlignment(c.outprefix + ".align3", final3, secrs, 3, a3Score, a1a2, c.linelimit);

    // Any het. InDel
    if (!bp.indelshift) {
      // Center on first SNP
      uint32_t reliableTracePos = findBestTraceSection(bc);
      bp.breakpoint = nearestSNP(c, bc, reliableTracePos);
    }

    // Variant Calling
    typedef std::vector<Variant> TVariants;
    TVariants var;
    if (c.callvariants) {
      now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Variant Calling" << std::endl;

      if (rs.forward) {
	callVariants(final1, allele1, var);
	callVariants(final2, allele2, var);
      } else {
	// Reverse complement
	std::string revPri(pri);
	reverseComplement(revPri);
	ReferenceSlice allele1Rev;
	_reverseReferenceSlize(allele1, allele1Rev);
	TAlign final1Rev;
	gotoh(revPri, allele1Rev.refslice, final1Rev, semiglobal, c.aliscore);
	callVariants(final1Rev, allele1Rev, var);
	std::string revSec(sec);
	reverseComplement(revSec);
	ReferenceSlice allele2Rev;
	_reverseReferenceSlize(allele2, allele2Rev);
	TAlign final2Rev;
	gotoh(revSec, allele2Rev.refslice, final2Rev, semiglobal, c.aliscore);
	callVariants(final2Rev, allele2Rev, var);
      }
	
      if ((c.annotatevariants) && (rs.filetype == 0)) {
	now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Variant Annotation (" << c.annotate << ")" << std::endl;
	std::string region = rs.chr + ":" + boost::lexical_cast<std::string>(rs.pos) + "-" + boost::lexical_cast<std::string>(rs.pos + rs.refslice.size());
	std::string response;
	if (!variantsInRegion(c, region, response)) {
	  std::vector<KnownVariation> kv;
	  int32_t numVar = parseKnownVariants(response, kv);
	  if (numVar > 0) {
	    annotateVariants(kv, var);
	  }
	} else {
	  std::cerr << "Warning: Variant annotation failed." << std::endl;
	  c.annotatevariants = false;
	}	  
      }

      // Sort variants
      std::sort(var.begin(), var.end());

      // VCF output
      vcfOutput(c, bc, var, rs);
    }
    
    // Json output
    traceAlleleAlignJsonOut(c, bc, tr, var, allele1, allele2, secrs, final1, final2, final3, dcp, a1Score, a2Score, a3Score, bp, a1a2);

    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Done." << std::endl;
    return 0;
  }

}

#endif
