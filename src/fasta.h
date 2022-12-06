#ifndef FASTA_H
#define FASTA_H

#include "abif.h"
#include "fmindex.h"

namespace tracy
{

#ifndef MAX_SINGLE_FASTA_SIZE
#define MAX_SINGLE_FASTA_SIZE 50000
#endif



  inline void
  _fixReferenceName(std::string& s) {
    // Disallow any weird characters
    boost::erase_all(s, "\\");
    boost::erase_all(s, ",");
    boost::erase_all(s, "'");
    boost::erase_all(s, "\"");
    boost::erase_all(s, "(");
    boost::erase_all(s, ")");
    boost::erase_all(s, "[");
    boost::erase_all(s, "]");
    boost::erase_all(s, "{");
    boost::erase_all(s, "}");
    boost::erase_all(s, "<");
    boost::erase_all(s, ">");
    boost::erase_all(s, ":");
    boost::erase_all(s, "\t");
    boost::erase_all(s, "\r");
    boost::erase_all(s, "#");
  }
  
  inline bool
  loadSingleFasta(std::string const& filename, std::string& faname, std::string& seq) {
    faname = "";
    std::string tmpfasta = "";
    std::ifstream fafile(filename.c_str());
    if (fafile.good()) {
      std::string line;
      while(std::getline(fafile, line)) {
	if (!line.empty()) {
	  if (line[0] == '>') {
	    if (!faname.empty()) {
	      std::cerr << "Only single-chromosome FASTA files are supported." << std::endl;
	      return false;
	    }
	    if (line.at(line.length() - 1) == '\r' ){
	      faname = line.substr(1, line.length() - 2);
	    } else {
	      faname = line.substr(1);
	    }
	  } else {
	    if (line.at(line.length() - 1) == '\r' ){
	      tmpfasta += boost::to_upper_copy(line.substr(0, line.length() - 1));
	    } else {
	      tmpfasta += boost::to_upper_copy(line);
	    }
	  }
	}
      }
      fafile.close();
    }
    // Check FASTA
    for(uint32_t k = 0; k < tmpfasta.size(); ++k)
      if ((tmpfasta[k] == 'A') || (tmpfasta[k] == 'C') || (tmpfasta[k] == 'G') || (tmpfasta[k] == 'T') || (tmpfasta[k] == 'N')) seq += tmpfasta[k];
    if (seq.size() != tmpfasta.size()) {
      std::cerr << "FASTA file contains nucleotides != [ACGTN]." << std::endl;
      return false;
    }

    // Fix FASTA sequence name for BCF output
    _fixReferenceName(faname);  // Replace special characters

    return true;
  }


  template<typename TConfig>
  inline void
  traceFastaOut(TConfig const& c, BaseCalls& bc, Trace const&) {
    // Output trace
    std::ofstream rfile(c.outfile.c_str());
    if (c.otype == "primary") {
      rfile << ">primary" << std::endl;
      for(uint32_t i = c.trimLeft; i < (bc.primary.size() - c.trimRight); ++i) rfile << bc.primary[i];
      rfile << std::endl;
    }
    else if (c.otype == "secondary") {
      rfile << ">secondary" << std::endl;
      for(uint32_t i = c.trimLeft; i < (bc.secondary.size() - c.trimRight); ++i) rfile << bc.secondary[i];
      rfile << std::endl;
    }
    else if (c.otype == "consensus") {
      rfile << ">consensus" << std::endl;
      for(uint32_t i = c.trimLeft; i < (bc.consensus.size() - c.trimRight); ++i) rfile << bc.consensus[i];
      rfile << std::endl;
    }
    rfile.close();  
  }

  template<typename TConfig>
  inline void
  traceFastqOut(TConfig const& c, BaseCalls& bc, Trace const& tr) {
    // Output trace
    std::ofstream rfile(c.outfile.c_str());
    if (c.otype == "primary") {
      rfile << "@primary" << std::endl;
      for(uint32_t i = c.trimLeft; i < (bc.primary.size() - c.trimRight); ++i) rfile << bc.primary[i];
      rfile << std::endl;
    }
    else if (c.otype == "secondary") {
      rfile << "@secondary" << std::endl;
      for(uint32_t i = c.trimLeft; i < (bc.secondary.size() - c.trimRight); ++i) rfile << bc.secondary[i];
      rfile << std::endl;
    }
    else if (c.otype == "consensus") {
      rfile << "@consensus" << std::endl;
      for(uint32_t i = c.trimLeft; i < (bc.consensus.size() - c.trimRight); ++i) rfile << bc.consensus[i];
      rfile << std::endl;
    }
    rfile << "+" << std::endl;
    
    typedef Trace::TValue TValue;
    uint32_t bcpos = 0;
    TValue idx = bc.bcPos[bcpos];
    for(int32_t i = 0; i < (int32_t) tr.traceACGT[0].size(); ++i) {
      if (idx == i) {
	if ((bcpos >= c.trimLeft) && (bcpos < (bc.primary.size() - c.trimRight))) rfile << (char) (bc.estQual[bcpos] + 33);
	if (bcpos < bc.bcPos.size() - 1) idx = bc.bcPos[++bcpos];
      }
    }
    rfile << std::endl;
    rfile.close();
  }

}

#endif
