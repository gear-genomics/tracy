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

#ifndef INDEX_H
#define INDEX_H

#include <sdsl/suffix_arrays.hpp>

#include <fstream>
#include <iomanip>

#include <boost/dynamic_bitset.hpp>
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

using namespace sdsl;

namespace tracy
{

  struct IndexConfig {
    boost::filesystem::path genome;
    boost::filesystem::path outfile;
  };

  int index(int argc, char** argv) {
    IndexConfig c;
  
    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("output,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("genome.fm9"), "output file")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.genome), "bgzipped genome file")
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
      std::cout << "Usage: tracy " << argv[0] << " [OPTIONS] genome.fa.gz" << std::endl;
      std::cout << visible_options << "\n";
      return -1;
    }
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "tracy ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    // Reference index
    csa_wt<> fm_index;
    
    // What kind of reference?
    std::ifstream ifile(c.genome.string().c_str(), std::ios::binary | std::ios::in);
    if (ifile.is_open()) {
      char fcode[4];
      ifile.seekg(0);
      ifile.read(fcode, 4);
      if (((uint8_t)fcode[0] == (uint8_t)0x1f) && ((uint8_t)fcode[1] == (uint8_t)0x8b)) {
	boost::filesystem::path dumpfile(c.outfile.string() + ".dump");
	std::string index_file = c.outfile.string();
	
	// Load FM index
	now = boost::posix_time::second_clock::local_time();
	if (!load_from_file(fm_index, index_file)) {
	  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Prepare FM-Index" << std::endl;
	  // Dump fasta
	  bool firstSeq = true;
	  std::ofstream tmpout(dumpfile.string().c_str());
	  std::ifstream file(c.genome.string().c_str(), std::ios_base::in | std::ios_base::binary);
	  boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
	  dataIn.push(boost::iostreams::gzip_decompressor());
	  dataIn.push(file);
	  std::istream instream(&dataIn);
	  std::string line;
	  while(std::getline(instream, line)) {
	    if ((!line.empty()) && (line[0] == '>')) {
	      if (!firstSeq) tmpout << std::endl;
	      else firstSeq = false;
	    } else {
	      tmpout << boost::to_upper_copy(line);
	    }
	  }
	  tmpout << std::endl;
	  file.close();
	  tmpout.close();
	  
	  now = boost::posix_time::second_clock::local_time();
	  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Create FM-Index" << std::endl;
	  
	  // Build index
	  construct(fm_index, dumpfile.string().c_str(), 1);
	  store_to_file(fm_index, index_file);
	  boost::filesystem::remove(dumpfile);
	}
      }
      ifile.close();
    }
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Done." << std::endl;
    
    return 0;
  }

}

#endif
