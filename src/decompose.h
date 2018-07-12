/*
============================================================================
Tracy: Trace File Handling
============================================================================
Copyright (C) 2017,2018 Tobias Rausch

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

#ifndef DECOMPOSE_H
#define DECOMPOSE_H

namespace tracy
{

  struct TraceBreakpoint {
    bool indelshift;
    bool traceleft;
    int32_t breakpoint;
    float bestDiff;
  };
  
  template<typename TConfig>
  inline void
  findBreakpoint(TConfig const& c, BaseCalls& bc, TraceBreakpoint& bp) {
    int32_t ncount = 0;
    for(uint32_t i = 0; ((i<c.kmer) && (i<bc.consensus.size())); ++i)
      if (bc.consensus[i] == 'N') ++ncount;
    std::vector<float> nratio;
    nratio.push_back((float)ncount / (float)c.kmer);  
    for(uint32_t i = c.kmer; i < bc.consensus.size(); ++i) {
      if (bc.consensus[i-c.kmer] == 'N') --ncount;
      if (bc.consensus[i] == 'N') ++ncount;
      nratio.push_back((float)ncount / (float)c.kmer);
    }
    float totalN = 0;
    for(uint32_t i = 0; i < nratio.size(); ++i) totalN += nratio[i];
    float leftSum = nratio[0];
    float rightSum = totalN - leftSum;
    bp.bestDiff = 0;
    bp.traceleft = true;
    bp.breakpoint = 0;
    for(uint32_t i = 1; i < nratio.size() - 1; ++i) {
      float right = rightSum / (float)(nratio.size() - i);
      float left = leftSum / (float)i;
      float diff = std::abs(right - left);
      if (diff > bp.bestDiff) {
	bp.breakpoint = i;
	bp.bestDiff = diff;
	if (left < right) bp.traceleft = true;
	else bp.traceleft = false;
      }
      leftSum += nratio[i];
      rightSum -= nratio[i];
    }
    bp.indelshift = true;
    // Forward breakpoint to first N
  for(uint32_t i = bp.breakpoint; i < bc.consensus.size(); ++i) {
    if (bc.consensus[i] == 'N') {
      bp.breakpoint = i;
      break;
    }
  }
  if ((bp.breakpoint <= c.trimLeft) || ((bc.consensus.size() - bp.breakpoint <= c.trimRight)) || (bp.bestDiff < 0.25)) {
    // No indel shift
    bp.indelshift = false;
    bp.breakpoint = bc.consensus.size() - c.trimRight - 1;
    bp.traceleft = true;
    bp.bestDiff = 0;
  } else {
    bp.indelshift = true;
  }    
}


  /*  

template<typename TConfig, typename TAlign>
inline void
plotAlignment(TConfig const& c, TAlign const& align, ReferenceSlice const& rs, int32_t key) {
  typedef typename TAlign::index TAIndex;
  uint32_t vi = 1;
  uint32_t ri = rs.pos;
  
  std::vector<char> alt;
  int32_t s = -1;
  int32_t e = -1;
  for(TAIndex j = 0; j< (TAIndex) align.shape()[1]; ++j) {
    if ((s == -1) && (align[1][j] != '-')) ++ri;
    if (align[0][j] != '-') {
      if (s == -1) s = j;
      e = j + 1;
      alt.push_back(align[0][j]);
    }
  }

  uint32_t riend = ri - 1;
  std::vector<char> ref;
  for(TAIndex j = s; j < (TAIndex) e; ++j) {
    if (align[1][j] != '-') {
      ref.push_back(align[1][j]);
      ++riend;
    }
  }


  boost::filesystem::path outalign(c.outprefix + ".align" + boost::lexical_cast<std::string>(key));
  std::ofstream ofile(outalign.string().c_str());
  uint32_t fald = c.linelimit + 14;
  ofile << ">Alt" << std::endl;
  for(uint32_t i = 0; i < alt.size(); ++i) {
    ofile << alt[i];
    if ((i+1) % fald == 0) ofile << std::endl;
  }
  if (alt.size() % fald != 0) ofile << std::endl;
  if (rs.forward) ofile << ">Ref " << rs.chr << ":" << ri << "-" << riend << " forward" << std::endl;
  else ofile << ">Ref " << rs.chr << ":" << rs.pos + rs.refslice.size() - (riend - rs.pos) + 1 << "-" << rs.pos + rs.refslice.size() - (ri - rs.pos) + 1 << " reversecomplement" << std::endl;
  for(uint32_t i = 0; i < ref.size(); ++i) {
    ofile << ref[i];
    if ((i+1)% fald == 0) ofile << std::endl;
  }
  if (ref.size() % fald != 0) ofile << std::endl;
  ofile << std::endl;
  ofile << "#";
  for(uint32_t i = 1; i < fald; ++i) ofile << "-";
  ofile << std::endl;
  ofile << std::endl;

  uint32_t blockcount = 0;
  while (s < e) {
    ofile << "Alt" << std::setw(10) << vi << ' ';
    for(TAIndex j = s; ((j < (TAIndex) e) && (j < s + c.linelimit)); ++j) {
      ofile << align[0][j];
      if (align[0][j] != '-') ++vi;
    }
    ofile << std::endl;
    ofile << "              ";
    for(TAIndex j = s; ((j < (TAIndex) e) && (j < s + c.linelimit)); ++j) {
      if (align[0][j] == align[1][j]) ofile << "|";
      else ofile << " ";
    }
    ofile << std::endl;
    if (rs.forward) ofile << "Ref" << std::setw(10) << ri << ' ';
    else ofile << "Ref" << std::setw(10) << rs.pos + rs.refslice.size() - (ri - rs.pos) + 1 << ' ';
    for(TAIndex j = s; ((j < (TAIndex) e) && (j < s + c.linelimit)); ++j) {
      ofile << align[1][j];
      if (align[1][j] != '-') ++ri;
    }
    ofile << std::endl;
    ofile << std::endl;
    s += c.linelimit;
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
 

  
template<typename TIterator, typename TValue>
inline void
getMedian(TIterator begin, TIterator end, TValue& median) {
  std::nth_element(begin, begin + (end - begin) / 2, end);
  median = *(begin + (end - begin) / 2);
}

template<typename TIterator, typename TValue>
inline void
getMAD(TIterator begin, TIterator end, TValue median, TValue& mad) {
  std::vector<TValue> absDev;
  for(;begin<end;++begin) 
    absDev.push_back(std::abs((TValue)*begin - median));
  getMedian(absDev.begin(), absDev.end(), mad);
}
 
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

inline char
phaseRefAllele(BaseCalls& bc, char const r, uint32_t varIndex) {
  if (bc.secondary[varIndex] == r) return bc.primary[varIndex];
  else if (bc.secondary[varIndex] == 'N') return 'N';
  else {
    if (bc.secondary[varIndex] == 'R') {
      if (r == 'A') return iupac(bc.primary[varIndex], 'G');
      else if (r == 'G') return iupac(bc.primary[varIndex], 'A');
    } else if (bc.secondary[varIndex] == 'Y') {
      if (r == 'C') return iupac(bc.primary[varIndex], 'T');
      else if (r == 'T') return iupac(bc.primary[varIndex], 'C');
    } else if (bc.secondary[varIndex] == 'S') {
      if (r == 'C') return iupac(bc.primary[varIndex], 'G');
      else if (r == 'G') return iupac(bc.primary[varIndex], 'C');
    } else if (bc.secondary[varIndex] == 'W') {
      if (r == 'A') return iupac(bc.primary[varIndex], 'T');
      else if (r == 'T') return iupac(bc.primary[varIndex], 'A');
    } else if (bc.secondary[varIndex] == 'K') {
      if (r == 'G') return iupac(bc.primary[varIndex], 'T');
      else if (r == 'T') return iupac(bc.primary[varIndex], 'G');
    } else if (bc.secondary[varIndex] == 'M') {
      if (r == 'A') return iupac(bc.primary[varIndex], 'C');
      else if (r == 'C') return iupac(bc.primary[varIndex], 'A');
    } else {
      return 'N';
    }
  }
  return 'N';
}


template<typename TConfig>
inline bool
findHomozygousBreakpoint(TConfig const& c, BaseCalls& bc, ReferenceSlice& rs) {
  // Homozygous mutation, estimate the breakpoint based on percent identity
  typedef boost::multi_array<char, 2> TAlign;
  TAlign align;
  AlignConfig<true, false> semiglobal;
  DnaScore<int> sc(5, -4, -10, -1);
  std::string consslice = trimmedCSeq(bc);
  gotoh(consslice, rs.refslice, align, semiglobal, sc);

  typedef typename TAlign::index TAIndex;
  TAIndex alignStart = 0;
  TAIndex alignEnd = 0;
  for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
    if ((align[0][j] != '-') && (align[1][j] != '-')) {
      alignStart = j;
      break;
    }
  }
  for(int32_t j = (int32_t) (align.shape()[1] - 1); j >= 0; --j) {
    if ((align[0][j] != '-') && (align[1][j] != '-')) {
      alignEnd = j;
      break;
    }
  }
  if (alignStart >= alignEnd) {
    std::cerr << "No valid alignment found between consensus and reference!" << std::endl;
    return false;
  }
  int32_t mismatch = 0;
  for(uint32_t j = alignStart; ((j<c.kmer) && (j<alignEnd)); ++j)
    if (align[0][j] != align[1][j]) ++mismatch;
  std::vector<float> mmratio;
  mmratio.push_back((float)mismatch/(float)c.kmer);
  for(uint32_t j = alignStart + c.kmer; j<alignEnd; ++j) {
    if (align[0][j-c.kmer] != align[1][j-c.kmer]) --mismatch;
    if (align[0][j] != align[1][j]) ++mismatch;
    mmratio.push_back((float)mismatch/(float)c.kmer);
  }
  float totalMM = 0;
  for(uint32_t i = 0; i < mmratio.size(); ++i) totalMM += mmratio[i];
  float leftSum = mmratio[0];
  float rightSum = totalMM - leftSum;
  float bestDiff = 0;
  bool traceleft = true;
  TAIndex bp = 0;
  for(uint32_t i = 1; i < mmratio.size() - 1; ++i) {
    float right = rightSum / (float)(mmratio.size() - i);
    float left = leftSum / (float)i;
    float diff = std::abs(right - left);
    if (diff > bestDiff) {
      bp = i;
      bestDiff = diff;
      if (left < right) traceleft = true;
      else traceleft = false;
    }
    leftSum += mmratio[i];
    rightSum -= mmratio[i];
  }
  // Find true consensus sequence breakpoint
  TAIndex varIndex = bc.ltrim;
  for(TAIndex j = alignStart; j < alignEnd; ++j) {
    if (j >= alignStart + bp) {
      // Forward breakpoint to first N or gap
      if ((align[0][j] == '-') || (align[0][j] == 'N')) break;
    }
    if (align[0][j] != '-') ++varIndex;
  }
  bc.breakpoint = varIndex;
  bc.indelshift = true;
  if (bestDiff < 0.25) {
    // Likely no hom. indel
    bc.indelshift = false;
    bc.breakpoint = bc.consensus.size() / 2;
    traceleft = true;
    bestDiff = 0;
  }

  //std::cout << bc.indelshift << ',' << bc.breakpoint << ',' << bc.consensus.size() << ',' << traceleft << std::endl;
  if (!traceleft) {
    bc.breakpoint = (uint16_t) (bc.consensus.size() - bc.breakpoint - 1);
    std::reverse(bc.consensus.begin(), bc.consensus.end());
    std::reverse(bc.primary.begin(), bc.primary.end());
    std::reverse(bc.secondary.begin(), bc.secondary.end());
    uint16_t tmptrim = bc.ltrim;
    bc.ltrim = bc.rtrim;
    bc.rtrim = tmptrim;
    for(uint32_t k = 0; k<4; ++k) {
      std::reverse(bc.peak[k].begin(), bc.peak[k].end());
      std::reverse(bc.pos[k].begin(), bc.pos[k].end());
    }
  }
  
  return true;
}  

 
template<typename TConfig>
inline bool
decomposeAlleles(TConfig const& c, BaseCalls& bc, ReferenceSlice& rs) {
  if (bc.ltrim >= bc.breakpoint) {
    std::cerr << "Breakpoint is inside the trimmed boundaries of the Sanger trace." << std::endl;
    return false;
  }
  // Align consensus to reference
  typedef boost::multi_array<char, 2> TAlign;
  TAlign align;
  AlignConfig<true, false> semiglobal;
  DnaScore<int> sc(5, -4, -10, -1);
  std::string consslice = bc.consensus.substr(bc.ltrim, (bc.breakpoint-bc.ltrim));
  gotoh(consslice, rs.refslice, align, semiglobal, sc);

  // Find breakpoint in alignment
  uint32_t varIndex = 0;
  uint32_t refPointer = 0;
  uint32_t alignIndex = 0;
  uint32_t vi = bc.ltrim;
  for(uint32_t j = 0; j < align.shape()[1]; ++j) {
    if (align[0][j] != '-') {
      if (align[1][j] != bc.primary[vi]) {
	char sec = phaseRefAllele(bc, align[1][j], vi);
	if (sec != 'N') {
	  bc.primary[vi] = align[1][j];
	  bc.secondary[vi] = sec;
	}
      }
      ++vi;
      if (vi == bc.breakpoint) {
	alignIndex = j;
	varIndex = vi;
	break;
      }
    }
    if (align[1][j] != '-') ++refPointer;
  }

  // Iterate possible deletion lengths
  std::vector<int32_t> fref;
  uint32_t maxdel = 2;
  if (rs.refslice.size() > (refPointer + bc.rtrim + 2)) maxdel = rs.refslice.size() - (refPointer + bc.rtrim);
  for(uint32_t del = 0; ((del < c.maxindel) && (del < maxdel / 2)); ++del) {
    int32_t failedref = 0;
    vi = varIndex;
    for(uint32_t j = alignIndex + del + 1; ((j < align.shape()[1]) && (vi < (bc.consensus.size() - bc.rtrim))); ++j, ++vi) {
      if (align[1][j] != bc.primary[vi]) {
	char sec = phaseRefAllele(bc, align[1][j], vi);
	if (sec == 'N') ++failedref;
      }
    }
    fref.push_back(failedref);
  }

  // Estimate cutoffs
  std::vector<int32_t> gm(fref);
  int32_t med = 0;
  getMedian(gm.begin(), gm.end(), med);
  int32_t mad = 0;
  getMAD(gm.begin(), gm.end(), med, mad);
  int32_t thres = 0;
  if (med > c.madc * mad) thres = med - c.madc * mad;
  if (thres < 10) thres = 10;
  //std::cout << thres << ',' << med << ','<< mad << std::endl;

  // Decompose using minimum deletion length
  std::vector<int32_t> deldecomp;
  for(uint32_t i = 0; i < fref.size(); ++i) {
    if (fref[i] < thres) {
      if ((i + 1 < fref.size()) && (2 * fref[i] < fref[i+1])) deldecomp.push_back(i);
      else if ((i > 0) && (2 * fref[i] < fref[i-1])) deldecomp.push_back(i);
    }
  }

  // Iterate possible insertion lengths
  std::vector<int32_t> fins;
  fins.push_back(fref[0]);
  uint32_t maxins = (int32_t) bc.consensus.size() - (int32_t) (bc.rtrim + bc.breakpoint);
  for(uint32_t ins = 1; ((ins < c.maxindel) && (ins < maxins / 2)); ++ins) {
    int32_t failedref = 0;
    vi = varIndex + ins;
    for(uint32_t j = alignIndex + 1; ((j < align.shape()[1]) && (vi < (bc.consensus.size() - bc.rtrim))); ++j, ++vi) {
      if (align[1][j] != bc.primary[vi]) {
	char sec = phaseRefAllele(bc, align[1][j], vi);
	if (sec == 'N') ++failedref;
      }
    }
    fins.push_back(failedref);
  }

  // Decompose using minimum insertion length
  std::vector<int32_t> insdecomp;
  for(uint32_t i = 0; i < fins.size(); ++i) {
    if (fins[i] < thres) {
      if ((i + 1 < fins.size()) && (2 * fins[i] < fins[i+1])) insdecomp.push_back(i);
      else if ((i > 0) && (2 * fins[i] < fins[i-1])) insdecomp.push_back(i);
    }
  }

  // Output decomposition table
  int32_t defins = 15;
  if ((deldecomp.empty()) && (insdecomp.empty())) defins = 50;
  for(uint32_t i = 0; i < insdecomp.size(); ++i)
    if (insdecomp[i] + 15 > defins) defins = insdecomp[i] + 15;
  if (defins > (int32_t) fins.size()) defins = fins.size();
  int32_t defdel = 15;
  if ((deldecomp.empty()) && (insdecomp.empty())) defdel = 50;
  for(uint32_t i = 0; i < deldecomp.size(); ++i)
    if (deldecomp[i] + 15 > defdel) defdel = deldecomp[i] + 15;
  if (defdel > (int32_t) fref.size()) defdel = fref.size();
  boost::filesystem::path outdecomp(c.outprefix + ".decomp");
  std::ofstream ofile(outdecomp.string().c_str());
  ofile << "indel\tdecomp" << std::endl;
  for(int32_t i = defdel - 1; i>=0; --i) ofile << (-1 * i) << "\t" << fref[i] << std::endl;
  for(int32_t i = 1; i<defins; ++i) ofile << (i) << "\t" << fins[i] << std::endl;
  ofile.close();

  // Actual decomposition
  if ((deldecomp.empty()) && (insdecomp.empty())) {
    // complex mutation
    int32_t bestIns = 0;
    int32_t bestDel = 0;
    int32_t bestFR = 1000;
    for(uint32_t ins = 0; ((ins < c.maxindel) && (ins < maxins / 2)); ++ins) {
      int32_t prevFailedRef = 0;
      for(uint32_t del = 0; ((del < c.maxindel) && (del <  maxdel / 2)); ++del) {
	int32_t failedref = 0;
	vi = varIndex + ins;
	for(uint32_t j = alignIndex + del + 1; ((j < align.shape()[1]) && (vi < (bc.consensus.size() - bc.rtrim))); ++j, ++vi) {
	  if (align[1][j] != bc.primary[vi]) {
	    char sec = phaseRefAllele(bc, align[1][j], vi);
	    if (sec == 'N') ++failedref;
	  }
	}
	if (2*failedref < prevFailedRef) {
	  if (failedref < bestFR) {
	    bestIns = ins;
	    bestDel = del;
	    bestFR = failedref;
	  }
	}
	prevFailedRef = failedref;
      }
    }
    if (bestFR != 1000) {
      std::cout << "Complex mutation, decomposition: ins: " << bestIns << ", del: " << bestDel << ", error: " << bestFR << std::endl;
      vi = varIndex + bestIns;
      for(uint32_t j = alignIndex + bestDel + 1; ((j < align.shape()[1]) && (vi < (bc.consensus.size() - bc.rtrim))); ++j, ++vi) {
	if (align[1][j] != bc.primary[vi]) {
	  char sec = phaseRefAllele(bc, align[1][j], vi);
	  if (sec != 'N') {
	    bc.primary[vi] = align[1][j];
	    bc.secondary[vi] = sec;
	  }
	}
      }
    } else {
      std::cout << "Allele decomposition failed, primary & secondary base calls unchanged." << std::endl;
    }
  } else {
    // take smallest deletion for decomposition
    if (!deldecomp.empty()) {
      std::sort(deldecomp.begin(), deldecomp.end());
      vi = varIndex;
      for(uint32_t j = alignIndex + deldecomp[0] + 1; ((j < align.shape()[1]) && (vi < (bc.consensus.size() - bc.rtrim))); ++j, ++vi) {
	if (align[1][j] != bc.primary[vi]) {
	  char sec = phaseRefAllele(bc, align[1][j], vi);
	  if (sec != 'N') {
	    bc.primary[vi] = align[1][j];
	    bc.secondary[vi] = sec;
	  }
	}
      }
    } else {
      std::sort(insdecomp.begin(), insdecomp.end());
      vi = varIndex + insdecomp[0];
      for(uint32_t j = alignIndex + 1; ((j < align.shape()[1]) && (vi < (bc.consensus.size() - bc.rtrim))); ++j, ++vi) {
	if (align[1][j] != bc.primary[vi]) {
	  char sec = phaseRefAllele(bc, align[1][j], vi);
	  if (sec != 'N') {
	    bc.primary[vi] = align[1][j];
	    bc.secondary[vi] = sec;
	  }
	}
      }
    }
  }
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
scanLeft(TFMIndex const& fm_index, std::string const& consensus, uint16_t const bestIdx, uint16_t const kmer, uint16_t const trim, THits& hits, bool unique) {
  int32_t ncount = 0;
  for(uint16_t i = bestIdx; ((i < bestIdx + kmer) && (i < consensus.size())); ++i)
      if (consensus[i] == 'N') ++ncount;
  for(int16_t k = bestIdx - 1; k>=trim; --k) {
    if (consensus[k+kmer] == 'N') --ncount;
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
	if (occs > 0) {
	  auto locations = locate(fm_index, seq.begin(), seq.end());
	  for(std::size_t m = 0; m < occs; ++m) {
	    hits.push_back(locations[m] - k);
	  }
	}
      }
    }
  }
}



template<typename TFMIndex, typename THits>
inline void
scanRight(TFMIndex const& fm_index, std::string const& consensus, uint16_t const bestIdx, uint16_t const kmer, uint16_t const trim, THits& hits, bool unique) {
  int32_t ncount = 0;
  for(uint16_t i = bestIdx; ((i < bestIdx + kmer) && (i < consensus.size())); ++i)
    if (consensus[i] == 'N') ++ncount;
  for(uint16_t k = bestIdx + kmer; k < (consensus.size() - trim); ++k) {
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
	if (occs > 0) {
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
  uint32_t minKmerSupport = 3;
  
  // Get sequence lengths
  std::vector<uint32_t> seqlen;
  faidx_t* fai = NULL;
  if (c.filetype) {
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
  scanLeft(fm_index, bc.consensus, bc.breakpoint, c.kmer, bc.ltrim, hitFwd, true);
  std::string rv = bc.consensus;
  reverseComplement(rv);
  scanRight(fm_index, rv, (uint16_t) (rv.size() - bc.breakpoint - 1), c.kmer, bc.rtrim, hitRev, true);
  
  // Select best orientation
  int64_t bestFwd;
  uint32_t freqFwd = findMaxFreq(hitFwd, bestFwd);
  int64_t bestRev;
  uint32_t freqRev = findMaxFreq(hitRev, bestRev);
  int64_t bestPos;
  if ((freqFwd >= minKmerSupport) && (freqFwd > 2*freqRev)) {
    rs.forward = true;
    rs.kmersupport = freqFwd;
    bestPos = bestFwd;
  } else if ((freqRev >= minKmerSupport) && (freqRev > 2*freqFwd)) {
    rs.forward = false;
    rs.kmersupport = freqRev;
    bestPos = bestRev;
  } else {
    // Try using non-unique matches
    hitFwd.clear();
    hitRev.clear();
    scanLeft(fm_index, bc.consensus, bc.breakpoint, c.kmer, bc.ltrim, hitFwd, false);
    scanRight(fm_index, rv, (uint16_t) (rv.size() - bc.breakpoint - 1), c.kmer, bc.rtrim, hitRev, false);
    freqFwd = findMaxFreq(hitFwd, bestFwd);
    freqRev = findMaxFreq(hitRev, bestRev);
    if ((freqFwd >= minKmerSupport) && (freqFwd > 2*freqRev)) {
      rs.forward = true;
      rs.kmersupport = freqFwd;
      bestPos = bestFwd;
    } else if ((freqRev >= minKmerSupport) && (freqRev > 2*freqFwd)) {
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
  if (!c.filetype) rs.chr = std::string(faidx_iseq(fai, refIndex));
  uint32_t chrpos = bestPos - cumsum;
  int32_t slen = -1;
  uint32_t slicestart = 0;
  uint32_t sliceend = seqlen[refIndex];
  if (bc.indelshift) {
    if (rs.forward) {
      if (chrpos > (uint32_t) (0.05 * (float) (bc.consensus.size()))) slicestart = chrpos - (int) (0.05 * (float) (bc.consensus.size()));
      uint32_t tmpend = chrpos + bc.consensus.size() + c.maxindel;
      if (tmpend < seqlen[refIndex]) sliceend = tmpend;
    } else {
      if (chrpos > c.maxindel) slicestart = chrpos - c.maxindel;
      uint32_t tmpend = chrpos + bc.consensus.size() + (int) (0.05 * (float) (bc.consensus.size()));
      if (tmpend < seqlen[refIndex]) sliceend = tmpend;
    }
  } else {
    // Homozygous mutation, search both sides
    if (chrpos > c.maxindel) slicestart = chrpos - c.maxindel;
    uint32_t tmpend = chrpos + bc.consensus.size() + c.maxindel;
    if (tmpend < seqlen[refIndex]) sliceend = tmpend;
  }
  if (!c.filetype) {
    rs.pos = slicestart;
    char* seq = faidx_fetch_seq(fai, rs.chr.c_str(), slicestart, sliceend, &slen);
    rs.refslice = boost::to_upper_copy(std::string(seq));
    if (seq != NULL) free(seq);
  }
  if (!rs.forward) reverseComplement(rs.refslice);
  //std::cout << rs.chr << "\t" << rs.pos << "\t" << rs.forward << std::endl;
  
  // Clean-up
  if (fai != NULL) fai_destroy(fai);
 
  return true;
}
  */
 
}

#endif
