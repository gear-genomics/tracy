#ifndef TEAL_H
#define TEAL_H

#include "abif.h"
#include "scf.h"
#include "trim.h"
#include "json.h"
#include "fasta.h"

namespace tracy {

  struct TealConfig {
    float pratio;
    uint16_t trimLeft;
    uint16_t trimRight;
    float trimStringency;
    std::string format;
    std::string otype;
    boost::filesystem::path tracein;
    boost::filesystem::path outfile;
  };

  int teal(int argc, char** argv) {
    TealConfig c;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("pratio,p", boost::program_options::value<float>(&c.pratio)->default_value(0.33), "peak ratio to call a base")
      ("format,f", boost::program_options::value<std::string>(&c.format)->default_value("json"), "output format [json|tsv|fasta|fastq]")
      ("otype,y", boost::program_options::value<std::string>(&c.otype)->default_value("primary"), "fasta/fastq sequence [primary|secondary|consensus]")
      ("trim,t", boost::program_options::value<float>(&c.trimStringency)->default_value(0), "trimming stringency [1:9], 0: use trimLeft and trimRight")
      ("trimLeft,q", boost::program_options::value<uint16_t>(&c.trimLeft)->default_value(0), "trim size left")
      ("trimRight,u", boost::program_options::value<uint16_t>(&c.trimRight)->default_value(0), "trim size right")
      ("output,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("out.json"), "basecalling output")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.tracein), "input trace file")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    
    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file"))) {
      std::cout << "Usage: tracy " << argv[0] << " [OPTIONS] trace.ab1" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "tracy ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;

    // Check input trace
    if (!(boost::filesystem::exists(c.tracein) && boost::filesystem::is_regular_file(c.tracein) && boost::filesystem::file_size(c.tracein))) {
      std::cerr << "Input trace file is missing: " << c.tracein.string() << std::endl;
      return 1;
    } else {
      // Read *.ab1 file
      Trace tr;
      int32_t ft = traceFormat(c.tracein.string());
      if (ft == 0) {
	if (!readab(c.tracein.string(), tr)) return -1;
      } else if (ft == 1) {
	if (!readscf(c.tracein.string(), tr)) return -1;
      } else {
	std::cerr << "Unknown trace file type!" << std::endl;
	return -1;
      }

      // Call bases
      BaseCalls bc;
      basecall(tr, bc, c.pratio);

      // Trace trimming
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

      // Write bases
      if (c.format == "tsv") traceTxtOut(c.outfile.string(), bc, tr, c.trimLeft, c.trimRight);
      else if (c.format == "fasta") traceFastaOut(c, bc, tr);
      else if (c.format == "fastq") traceFastqOut(c, bc, tr);
      else traceJsonOut(c.outfile.string(), bc, tr);
    }

    // End
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
    return 0;
  }

}

#endif
