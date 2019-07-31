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

#ifndef FMINDEX_H
#define FMINDEX_H

#include <boost/progress.hpp>

#include <htslib/faidx.h>

#include "fasta.h"

namespace tracy
{

  inline void
  reverseComplement(std::string& sequence) {
    std::string rev = boost::to_upper_copy(std::string(sequence.rbegin(), sequence.rend()));
    std::size_t i = 0;
    for(std::string::iterator revIt = rev.begin(); revIt != rev.end(); ++revIt, ++i) {
      switch (*revIt) {
      case 'A': sequence[i]='T'; break;
      case 'C': sequence[i]='G'; break;
      case 'G': sequence[i]='C'; break;
      case 'T': sequence[i]='A'; break;
      case 'N': sequence[i]='N'; break;
      default: break;
      }
    }
  }
  

  struct ReferenceSlice {
    bool forward;
    int32_t filetype;   //-1: failure, 0: *fa.gz, 1: *.fa, 2: *.ab1
    uint32_t kmersupport;
    uint32_t pos;
    std::string chr;
    std::string refslice;

    ReferenceSlice() : forward(true), filetype(-1), kmersupport(0), pos(0), chr(""), refslice("") {}
  };


  inline void
  _reverseReferenceSlize(ReferenceSlice const& in, ReferenceSlice& out) {
    out.forward = !in.forward;
    out.filetype = in.filetype;
    out.kmersupport = in.kmersupport;
    out.pos = in.pos;
    out.chr = in.chr;
    out.refslice = in.refslice;
    reverseComplement(out.refslice);
  }

  struct TraceBreakpoint {
    bool indelshift;
    bool traceleft;
    uint32_t breakpoint;
    float bestDiff;
  };


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

  inline int32_t     // -1: failure, 0: Indexed genome, 1: fasta file, 2: trace
  genomeType(std::string const& path) {
    std::ifstream ifile(path.c_str(), std::ios::binary | std::ios::in);
    if (ifile.is_open()) {
      char fcode[4];
      ifile.seekg(0);
      ifile.read(fcode, 4);
      ifile.close();
      if (((uint8_t)fcode[0] == (uint8_t)0x1f) && ((uint8_t)fcode[1] == (uint8_t)0x8b)) return 0; // Gzipped fasta, big reference genome
      else if (traceFormat(path) >= 0) return 2; // Trace file
      else if (fcode[0] == '>') return 1; // Single-fasta file
    }
    return -1;
  }
    
  
  template<typename TConfig, typename TFMIdx>
  inline bool
  loadFMIdx(TConfig const& c, ReferenceSlice& rs, TFMIdx& fm_index) {
    
    // What kind of reference?
    std::ifstream ifile(c.genome.c_str(), std::ios::binary | std::ios::in);
    if (ifile.is_open()) {
      char fcode[4];
      ifile.seekg(0);
      ifile.read(fcode, 4);
      if (((uint8_t)fcode[0] == (uint8_t)0x1f) && ((uint8_t)fcode[1] == (uint8_t)0x8b)) {
	// Gzipped fasta
	rs.filetype = 0;
	boost::filesystem::path op = c.genome.parent_path() / c.genome.stem();
	boost::filesystem::path outfile(op.string() + ".dump");
	std::string index_file = op.string() + ".fm9";

	// Load FM index
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Load FM-Index" << std::endl;

	// Check if that's an old index file
	if (!load_from_checked_file(fm_index, index_file)) {
	  std::cerr << "Old index data format. Please rebuild your reference genome index!" << std::endl;
	  return -1;
	}
      } else if (traceFormat(c.genome.string()) >= 0) {
	rs.filetype = 2;
	
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Load ab1 wildtype" << std::endl;
	Trace wt;
	int32_t ft = traceFormat(c.genome.string());
	if (ft == 0) {
	  if (!readab(c.genome.string(), wt)) return -1;
	} else if (ft == 1) {
	  if (!readscf(c.genome.string(), wt)) return -1;
	} else {
	  std::cerr << "Unknown trace file type!" << std::endl;
	  return -1;
	}
	BaseCalls wtbc;
	basecall(wt, wtbc, c.pratio);
	rs.chr = "wildtype";
	rs.refslice = wtbc.primary;
	construct_im(fm_index, rs.refslice.c_str(), 1);
      }
      else if (fcode[0] == '>') {
	rs.filetype = 1;
	
	// Single FASTA file
	rs.chr = "";
	std::string tmpfasta = "";
	boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
	std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Load FASTA reference" << std::endl;
	std::ifstream fafile(c.genome.string().c_str());
	if (fafile.good()) {
	  std::string line;
	  while(std::getline(fafile, line)) {
	    if (!line.empty()) {
	      if (line[0] == '>') {
		if (!rs.chr.empty()) {
		  std::cerr << "Only single-chromosome FASTA files are supported. If you have a multi-FASTA file please use bgzip and index the FASTA file with samtools faidx!" << std::endl;
		  return false;
		}
                if (line.at(line.length() - 1) == '\r' ){
                  rs.chr = line.substr(1, line.length() - 2);
                } else {
		  rs.chr = line.substr(1);
                }
		_fixReferenceName(rs.chr);  // Replace special characters
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
	rs.refslice = "";
	for(uint32_t k = 0; k < tmpfasta.size(); ++k)
	  if ((tmpfasta[k] == 'A') || (tmpfasta[k] == 'C') || (tmpfasta[k] == 'G') || (tmpfasta[k] == 'T') || (tmpfasta[k] == 'N')) rs.refslice += tmpfasta[k];
	if (rs.refslice.size() != tmpfasta.size()) {
	  std::cerr << "FASTA file contains nucleotides != [ACGTN]." << std::endl;
	  return false;
	}
	construct_im(fm_index, rs.refslice.c_str(), 1);
      } else {
	std::cerr << "Couldn't recognize reference file format!" << std::endl;
	return false;
      }
    }
    ifile.close();
    return true;
  }


  template<typename THits, typename TValue>
  inline uint32_t
  findMaxFreq(THits& hits, TValue& gpos) {
    if (hits.empty()) {
      gpos = 0;
      return 0;
    }
    std::sort(hits.begin(), hits.end());
    TValue ihit = hits[0];
    int32_t freq = 1;
    int32_t bestFreq = 1;
    gpos = ihit;
    for(uint32_t i = 1; i<hits.size(); ++i) {
      if (hits[i] == ihit) {
	++freq;
	if (freq > bestFreq) {
	  gpos = ihit;
	  bestFreq = freq;
	}
      } else {
	ihit = hits[i];
	freq = 1;
      }
    }
    return bestFreq;
  }
  

  
  
  template<typename TFMIndex, typename THits>
  inline void
  scanSequence(TFMIndex const& fm_index, std::string const& consensus, uint16_t const trimLeft, uint16_t const trimRight, uint16_t const kmer, THits& hits, bool unique) {
    int32_t ncount = 0;
    for(uint16_t i = trimLeft; ((i < trimLeft + kmer) && (i < consensus.size())); ++i)
      if (consensus[i] == 'N') ++ncount;
    for(uint16_t k = trimLeft + kmer; ((k < (consensus.size() - trimRight)) && (k < consensus.size())); ++k) {
      if (consensus[k-kmer] == 'N') --ncount;
      if (consensus[k] == 'N') ++ncount;
      if (ncount == 0) {
	std::string seq = consensus.substr(k, kmer);
	std::size_t occs = sdsl::count(fm_index, seq.begin(), seq.end());
	if (unique) {
	  if (occs == 1) {
	    auto locations = locate(fm_index, seq.begin(), seq.end());
	    hits.push_back(locations[0] - k);
	  }
	} else {
	  // Exclude ubiquitously mapping k-mers
	  if ((occs > 0) && (occs < 1000)) {
	    auto locations = locate(fm_index, seq.begin(), seq.end());
	    for(std::size_t m = 0; m < occs; ++m) {
	      hits.push_back(locations[m] - k);
	    }
	  }
	}
      }
    }
  }


  
  template<typename TConfig, typename TFMIndex>
  inline bool
  getReferenceSlice(TConfig const& c, TFMIndex const& fm_index, BaseCalls const& bc, ReferenceSlice& rs) {
  
    // Get sequence lengths
    std::vector<uint32_t> seqlen;
    faidx_t* fai = NULL;
    if (rs.filetype) {
      seqlen.push_back(rs.refslice.size());
    } else {
      fai = fai_load(c.genome.string().c_str());
      seqlen.resize(faidx_nseq(fai));
      for(int32_t refIndex = 0; refIndex < faidx_nseq(fai); ++refIndex) {
	std::string seqname(faidx_iseq(fai, refIndex));
	seqlen[refIndex] = faidx_seq_len(fai, seqname.c_str()) + 1;
      }
    }

    // Fwd and rev index search
    std::vector<int64_t> hitFwd;
    std::vector<int64_t> hitRev;
    scanSequence(fm_index, bc.consensus, c.trimLeft, c.trimRight, c.kmer, hitFwd, true);
    std::string rv = bc.consensus;
    reverseComplement(rv);
    scanSequence(fm_index, rv, c.trimRight, c.trimLeft, c.kmer, hitRev, true);
  
    // Select best orientation
    int64_t bestFwd;
    uint32_t freqFwd = findMaxFreq(hitFwd, bestFwd);
    int64_t bestRev;
    uint32_t freqRev = findMaxFreq(hitRev, bestRev);
    int64_t bestPos;
    if ((freqFwd >= c.minKmerSupport) && (freqFwd > 2*freqRev)) {
      rs.forward = true;
      rs.kmersupport = freqFwd;
      bestPos = bestFwd;
    } else if ((freqRev >= c.minKmerSupport) && (freqRev > 2*freqFwd)) {
      rs.forward = false;
      rs.kmersupport = freqRev;
      bestPos = bestRev;
    } else {
      // Try using non-unique matches
      hitFwd.clear();
      hitRev.clear();
      scanSequence(fm_index, bc.consensus, c.trimLeft, c.trimRight, c.kmer, hitFwd, false);
      scanSequence(fm_index, rv, c.trimRight, c.trimLeft, c.kmer, hitRev, false);
      freqFwd = findMaxFreq(hitFwd, bestFwd);
      freqRev = findMaxFreq(hitRev, bestRev);
      if ((freqFwd >= c.minKmerSupport) && (freqFwd > 2*freqRev)) {
	rs.forward = true;
	rs.kmersupport = freqFwd;
	bestPos = bestFwd;
      } else if ((freqRev >= c.minKmerSupport) && (freqRev > 2*freqFwd)) {
	rs.forward = false;
	rs.kmersupport = freqRev;
	bestPos = bestRev;
      } else {
	std::cerr << "Couldn't anchor the Sanger trace in the selected reference genome." << std::endl;
	return false;
      }
    }

    // Get initial ref slice
    int64_t cumsum = 0;
    uint32_t refIndex = 0;
    for(; bestPos >= cumsum + seqlen[refIndex]; ++refIndex) cumsum += seqlen[refIndex];
    if (!rs.filetype) rs.chr = std::string(faidx_iseq(fai, refIndex));
    uint32_t chrpos = bestPos - cumsum;
    int32_t slen = -1;
    uint32_t slicestart = 0;
    uint32_t sliceend = seqlen[refIndex];
    if (chrpos > c.maxindel) slicestart = chrpos - c.maxindel;
    uint32_t tmpend = chrpos + bc.consensus.size() + c.maxindel;
    if (tmpend < seqlen[refIndex]) sliceend = tmpend;
    if (!rs.filetype) {
      rs.pos = slicestart;
      char* seq = faidx_fetch_seq(fai, rs.chr.c_str(), slicestart, sliceend, &slen);
      rs.refslice = boost::to_upper_copy(std::string(seq));
      if (seq != NULL) free(seq);
    }
    if (!rs.forward) reverseComplement(rs.refslice);
    //std::cerr << rs.chr << "\t" << rs.pos << "\t" << rs.forward << std::endl;
    
    // Clean-up
    if (fai != NULL) fai_destroy(fai);
    
    return true;
  }


  template<typename TAlign>
  inline void
  plotAlignment(std::string const& filename, TAlign const& align, ReferenceSlice const& rs, int32_t const key, int32_t const score, std::pair<double, double> const& a1a2, uint32_t const linelimit) {
    typedef typename TAlign::index TAIndex;
    int32_t ri = rs.pos + 1;
    int32_t riend = rs.pos + rs.refslice.size();
    int32_t vi = 1;
    
    uint32_t fald = linelimit + 14;
    std::ofstream ofile(filename.c_str());
    if (key == 0) ofile << ">Alt" << std::endl;
    else if (key == 1) ofile << ">Alt1 (Estimated allelic Fraction: " << a1a2.first << ")" << std::endl;
    else if (key == 2) ofile << ">Alt2 (Estimated allelic Fraction: " << a1a2.second << ")" << std::endl;
    else ofile << ">Alt1 (Estimated allelic Fraction: " << a1a2.first << ")" << std::endl;
    int32_t count = 0;
    for(TAIndex j = 0; j< (TAIndex) align.shape()[1]; ++j) {
      if (align[0][j] != '-')  {
	ofile << align[0][j];
	if ((count+1) % fald == 0) ofile << std::endl;
	++count;
      }
    }
    if (count % fald != 0) ofile << std::endl;
    if (key != 3) {
      if (rs.forward) ofile << ">Ref " << rs.chr << ":" << ri << "-" << riend << " forward" << std::endl;
      else ofile << ">Ref " << rs.chr << ":" << rs.pos + rs.refslice.size() - (riend - rs.pos) + 1 << "-" << rs.pos + rs.refslice.size() - (ri - rs.pos) + 1 << " reversecomplement" << std::endl;
    } else {
      ofile << ">Alt2 (Estimated allelic Fraction: " << a1a2.second << ")" << std::endl;
    }
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
    
    uint32_t blockcount = 0;
    int32_t s = 0;
    int32_t e = (TAIndex) align.shape()[1];
    while (s < e) {
      if (key != 3) ofile << "Alt" << std::setw(10) << vi << ' ';
      else ofile << "Alt1" << std::setw(9) << vi << ' ';
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
      if (key != 3) {
	if (rs.forward) ofile << "Ref" << std::setw(10) << ri << ' ';
	else ofile << "Ref" << std::setw(10) << rs.pos + rs.refslice.size() - (ri - rs.pos) + 1 << ' ';
      } else {
	ofile << "Alt2" << std::setw(9) << ri << ' ';
      }
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

  template<typename TAlign>
  inline void
  plotAlignment(std::string const& filename, TAlign const& align, ReferenceSlice const& rs, int32_t const score, uint32_t const linelimit) {
    plotAlignment(filename, align, rs, 0, score, std::make_pair(0, 0), linelimit);
  }

  template<typename TConfig, typename TAlign>
  inline void
  trimReferenceSlice(TConfig const& c, TAlign const& align, ReferenceSlice& rs) {
    typedef typename TAlign::index TAIndex;
    uint32_t ri = 0;
    int32_t s = -1;
    int32_t e = -1;
    for(TAIndex j = 0; j< (TAIndex) align.shape()[1]; ++j) {
      if (align[0][j] != '-') {
	if (s == -1) s = j;
	e = j + 1;
      }
      if ((s == -1) && (align[1][j] != '-')) ++ri;
    }
    uint32_t risize = 0;
    for(TAIndex j = s; j < (TAIndex) e; ++j) {
      if (align[1][j] != '-') ++risize;
    }
    if (ri >= c.trimLeft) {
      ri -= c.trimLeft;
      risize += c.trimLeft;
    }
    if (ri + risize + c.trimRight < rs.refslice.size()) risize += c.trimRight;
    int32_t oldlen = rs.refslice.size();
    rs.refslice = rs.refslice.substr(ri, risize);
    if (rs.forward) rs.pos += ri;
    else {
      int32_t offset = oldlen - (int32_t) ri - (int32_t) risize;
      if (offset < 0) {
	std::cerr << "Warning: Offset smaller than zero!" << std::endl;
      } else {
	rs.pos += offset;
      }
    }
  }
  


}

#endif
