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

  template<typename TProfile>
  inline void
  findBreakpoint(TProfile const& ptrace, TraceBreakpoint& bp) {
    // Compute signal vector
    std::vector<double> sigratio;
    for(uint32_t j = 0; j<ptrace.shape()[1]; ++j) {
      double best = 0.001;
      double sndBest = 0.001;
      for(uint32_t i = 0; i<ptrace.shape()[0]; ++i) {
	if (ptrace[i][j] > best) {
	  sndBest = best;
	  best = ptrace[i][j];
	} else if (ptrace[i][j] > sndBest) {
	  sndBest = ptrace[i][j];
	}
      }
      sigratio.push_back(best - sndBest);
    }

    // Find best breakpoint
    bp.bestDiff = 0;
    bp.traceleft = true;
    bp.breakpoint = 0;
    uint32_t minWindow = 25;
    for(uint32_t i = minWindow; i < sigratio.size() - minWindow; ++i) {
      double leftSum = 0;
      for(uint32_t k = i-minWindow; k < i; ++k) leftSum += sigratio[k];
      double left = leftSum / (double) minWindow;
      double rightSum = 0;
      for(uint32_t k = i; k < i + minWindow; ++k) rightSum += sigratio[k];
      double right = rightSum / (double) minWindow;
      double diff = std::abs(right - left);
      if (diff > bp.bestDiff) {
	bp.breakpoint = i;
	bp.bestDiff = diff;
	if (left < right) bp.traceleft = false;
	else bp.traceleft = true;
      }
    }
    bp.indelshift = true;
    if (bp.bestDiff < 0.25) {
      // No indel shift
      bp.indelshift = false;
      bp.breakpoint = ptrace.shape()[1];
      bp.traceleft = true;
      bp.bestDiff = 0;
    }
  }


  template<typename TAlign>
  inline bool
  findHomozygousBreakpoint(TAlign& align, TraceBreakpoint& bp) {
    // Homozygous mutation, estimate the breakpoint based on percent identity
    typedef typename TAlign::index TAIndex;
    TAIndex alignStart = 0;
    TAIndex alignEnd = 0;
    TAIndex varIndex = 0;
    for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
      if ((align[0][j] != '-') && (align[1][j] != '-')) {
	alignStart = j;
	break;
      }
      if (align[0][j] != '-') ++varIndex;
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


    // Find best breakpoint
    bp.bestDiff = 0;
    bp.traceleft = true;
    bp.breakpoint = 0;
    uint32_t minWindow = 25;
    for(uint32_t i = alignStart; i < alignStart + minWindow; ++i) {
      if (align[0][i] != '-') ++varIndex;
    }
    for(uint32_t i = alignStart + minWindow; i < alignEnd - minWindow; ++i) {
      if (align[0][i] != '-') ++varIndex;
      double leftSum = 0;
      for(uint32_t k = i-minWindow; k < i; ++k) {
	if (align[0][k] != align[1][k]) ++leftSum;
      }
      double left = leftSum / (double) minWindow;
      double rightSum = 0;
      for(uint32_t k = i; k < i + minWindow; ++k) {
	if (align[0][k] != align[1][k]) ++rightSum;
      }
      double right = rightSum / (double) minWindow;
      double diff = std::abs(right - left);
      if (diff > bp.bestDiff) {
	bp.breakpoint = varIndex;
	bp.bestDiff = diff;
	if (left < right) bp.traceleft = true;
	else bp.traceleft = false;
      }
      //std::cerr << varIndex << ':' << diff << std::endl;
    }
    bp.indelshift = true;
    if (bp.bestDiff < 0.25) {
      // No indel shift
      bp.indelshift = false;
      bp.breakpoint = varIndex;
      bp.traceleft = true;
      bp.bestDiff = 0;
    }
    return true;
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


 
  template<typename TConfig, typename TAlign, typename TDecomp>
  inline bool
  decomposeAlleles(TConfig const& c, TAlign const& align, BaseCalls& bc, TraceBreakpoint bp, ReferenceSlice& rs, TDecomp& dcp) {
    int32_t ltrim = c.trimLeft;
    int32_t rtrim = c.trimRight;
    
    // Find breakpoint in alignment
    uint32_t varIndex = 0;
    uint32_t refPointer = 0;
    uint32_t alignIndex = 0;
    uint32_t vi = ltrim;
    bp.breakpoint += ltrim;
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
	if (vi == bp.breakpoint) {
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
    if (rs.refslice.size() > (refPointer + rtrim + 2)) maxdel = rs.refslice.size() - (refPointer + rtrim);
    for(uint32_t del = 0; ((del < c.maxindel) && (del < maxdel / 2)); ++del) {
      int32_t failedref = 0;
      vi = varIndex;
      for(uint32_t j = alignIndex + del + 1; ((j < align.shape()[1]) && (vi < (bc.consensus.size() - rtrim))); ++j, ++vi) {
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
    uint32_t maxins = (int32_t) bc.consensus.size() - (int32_t) (rtrim + bp.breakpoint);
    for(uint32_t ins = 1; ((ins < c.maxindel) && (ins < maxins / 2)); ++ins) {
      int32_t failedref = 0;
      vi = varIndex + ins;
      for(uint32_t j = alignIndex + 1; ((j < align.shape()[1]) && (vi < (bc.consensus.size() - rtrim))); ++j, ++vi) {
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
    for(int32_t i = defdel - 1; i>=0; --i) dcp.push_back(std::make_pair((-1 * i), fref[i]));
    for(int32_t i = 1; i<defins; ++i) dcp.push_back(std::make_pair(i, fins[i]));
    
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
	  for(uint32_t j = alignIndex + del + 1; ((j < align.shape()[1]) && (vi < (bc.consensus.size() - rtrim))); ++j, ++vi) {
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
	for(uint32_t j = alignIndex + bestDel + 1; ((j < align.shape()[1]) && (vi < (bc.consensus.size() - rtrim))); ++j, ++vi) {
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
	for(uint32_t j = alignIndex + deldecomp[0] + 1; ((j < align.shape()[1]) && (vi < (bc.consensus.size() - rtrim))); ++j, ++vi) {
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
	for(uint32_t j = alignIndex + 1; ((j < align.shape()[1]) && (vi < (bc.consensus.size() - rtrim))); ++j, ++vi) {
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

  inline void
  generateSecondaryDecomposed(Trace const& tr, BaseCalls& bc) {
    bc.secDecompose.resize(bc.secondary.size());
    for(uint32_t i = 0; ((i < bc.primary.size()) && (i < bc.secondary.size())); ++i) {
      if (bc.primary[i] == bc.secondary[i]) bc.secDecompose[i] = bc.primary[i];
      else {
	if (!isAmbiguous(bc.secondary[i])) bc.secDecompose[i] = bc.secondary[i];
	else {
	  uint32_t tracePos = bc.bcPos[i];
	  // Select higher peak
	  if (bc.secondary[i] == 'R') {
	    if (tr.traceACGT[0][tracePos] > tr.traceACGT[2][tracePos]) bc.secDecompose[i] = 'A';
	    else bc.secDecompose[i] = 'G';
	  } else if (bc.secondary[i] == 'Y') {
	    if (tr.traceACGT[1][tracePos] > tr.traceACGT[3][tracePos]) bc.secDecompose[i] = 'C';
	    else bc.secDecompose[i] = 'T';
	  } else if (bc.secondary[i] == 'S') {
	    if (tr.traceACGT[1][tracePos] > tr.traceACGT[2][tracePos]) bc.secDecompose[i] = 'C';
	    else bc.secDecompose[i] = 'G';
	  } else if (bc.secondary[i] == 'W') {
	    if (tr.traceACGT[0][tracePos] > tr.traceACGT[3][tracePos]) bc.secDecompose[i] = 'A';
	    else bc.secDecompose[i] = 'T';
	  } else if (bc.secondary[i] == 'K') {
	    if (tr.traceACGT[2][tracePos] > tr.traceACGT[3][tracePos]) bc.secDecompose[i] = 'G';
	    else bc.secDecompose[i] = 'T';
	  } else if (bc.secondary[i] == 'M') {
	    if (tr.traceACGT[0][tracePos] > tr.traceACGT[1][tracePos]) bc.secDecompose[i] = 'A';
	    else bc.secDecompose[i] = 'C';
	  } else bc.secDecompose[i] = 'N';
	}
      }
    }
  }
  
  template<typename TConfig>
  inline std::pair<double, double>
  allelicFraction(TConfig const& c, Trace const& tr, BaseCalls const& bc) {
    std::string pri = trimmedSeq(bc.primary, c.trimLeft, c.trimRight);
    std::string sec = trimmedSeq(bc.secDecompose, c.trimLeft, c.trimRight);
    uint32_t diffnuc = 0;
    for(uint32_t i = 0; ((i < pri.size()) && (i < sec.size())); ++i) {
      if (pri[i] != sec[i]) ++diffnuc;
    }

    double bestSSE = 0;
    double bestI = 0.5;
    double bestJ = 0.5;
    double bestK = 0;
    double bestL = 0;
    if (diffnuc) {
      // Get the max. 4 possible alleles
      uint32_t nucpos = 0;
      typedef boost::multi_array<double, 2> TProfile;
      TProfile tp(boost::extents[4][diffnuc]);  // A,C,G,T
      TProfile prip(boost::extents[4][diffnuc]);  // A,C,G,T
      TProfile secp(boost::extents[4][diffnuc]);  // A,C,G,T
      TProfile terp(boost::extents[4][diffnuc]);  // A,C,G,T
      TProfile quap(boost::extents[4][diffnuc]);  // A,C,G,T
      for(uint32_t i = 0; i<4; ++i) {
	for(uint32_t j = 0; j<diffnuc; ++j) {
	  tp[i][j] = 0;
	  prip[i][j] = 0;
	  secp[i][j] = 0;
	  terp[i][j] = 0;
	  quap[i][j] = 0;
	}
      }
      for(uint32_t i = 0; ((i < pri.size()) && (i < sec.size())); ++i) {
	if (pri[i] != sec[i]) {
	  uint32_t tracePos = bc.bcPos[i+c.trimLeft];
	  double sigsum = tr.traceACGT[0][tracePos] + tr.traceACGT[1][tracePos] + tr.traceACGT[2][tracePos] + tr.traceACGT[3][tracePos];
	  for(uint32_t k = 0; k<4; ++k) tp[k][nucpos] = (double) (tr.traceACGT[k][tracePos]) / sigsum;	  
	  //std::cout << pri[i] << ',' << sec[i] << ':' << tr.traceACGT[0][tracePos] << ',' << tr.traceACGT[1][tracePos] << ',' << tr.traceACGT[2][tracePos] << ',' << tr.traceACGT[3][tracePos] << std::endl;
	  if ((pri[i] == 'A') && (sec[i] =='C')) {
	    prip[0][nucpos] = 1;
	    secp[1][nucpos] = 1;
	    if (tr.traceACGT[2][tracePos] > tr.traceACGT[3][tracePos]) {
	      terp[2][nucpos] = 1;
	      quap[3][nucpos] = 1;
	    } else {
	      terp[3][nucpos] = 1;
	      quap[2][nucpos] = 1;
	    }
	  } else if ((pri[i] == 'A') && (sec[i] =='G')) {
	    prip[0][nucpos] = 1;
	    secp[2][nucpos] = 1;
	    if (tr.traceACGT[1][tracePos] > tr.traceACGT[3][tracePos]) {
	      terp[1][nucpos] = 1;
	      quap[3][nucpos] = 1;
	    } else {
	      terp[3][nucpos] = 1;
	      quap[1][nucpos] = 1;
	    }
	  } else if ((pri[i] == 'A') && (sec[i] =='T')) {
	    prip[0][nucpos] = 1;
	    secp[3][nucpos] = 1;
	    if (tr.traceACGT[1][tracePos] > tr.traceACGT[2][tracePos]) {
	      terp[1][nucpos] = 1;
	      quap[2][nucpos] = 1;
	    } else {
	      terp[2][nucpos] = 1;
	      quap[1][nucpos] = 1;
	    }
	  } else if ((pri[i] == 'C') && (sec[i] =='A')) {
	    prip[1][nucpos] = 1;
	    secp[0][nucpos] = 1;
	    if (tr.traceACGT[2][tracePos] > tr.traceACGT[3][tracePos]) {
	      terp[2][nucpos] = 1;
	      quap[3][nucpos] = 1;
	    } else {
	      terp[3][nucpos] = 1;
	      quap[2][nucpos] = 1;
	    }
	  } else if ((pri[i] == 'C') && (sec[i] =='G')) {
	    prip[1][nucpos] = 1;
	    secp[2][nucpos] = 1;
	    if (tr.traceACGT[0][tracePos] > tr.traceACGT[3][tracePos]) {
	      terp[0][nucpos] = 1;
	      quap[3][nucpos] = 1;
	    } else {
	      terp[3][nucpos] = 1;
	      quap[0][nucpos] = 1;
	    }
	  } else if ((pri[i] == 'C') && (sec[i] =='T')) {
	    prip[1][nucpos] = 1;
	    secp[3][nucpos] = 1;
	    if (tr.traceACGT[0][tracePos] > tr.traceACGT[2][tracePos]) {
	      terp[0][nucpos] = 1;
	      quap[2][nucpos] = 1;
	    } else {
	      terp[2][nucpos] = 1;
	      quap[0][nucpos] = 1;
	    }
	  } else if ((pri[i] == 'G') && (sec[i] =='A')) {
	    prip[2][nucpos] = 1;
	    secp[0][nucpos] = 1;
	    if (tr.traceACGT[1][tracePos] > tr.traceACGT[3][tracePos]) {
	      terp[1][nucpos] = 1;
	      quap[3][nucpos] = 1;
	    } else {
	      terp[3][nucpos] = 1;
	      quap[1][nucpos] = 1;
	    }
	  } else if ((pri[i] == 'G') && (sec[i] =='C')) {
	    prip[2][nucpos] = 1;
	    secp[1][nucpos] = 1;
	    if (tr.traceACGT[0][tracePos] > tr.traceACGT[3][tracePos]) {
	      terp[0][nucpos] = 1;
	      quap[3][nucpos] = 1;
	    } else {
	      terp[3][nucpos] = 1;
	      quap[0][nucpos] = 1;
	    }
	  } else if ((pri[i] == 'G') && (sec[i] =='T')) {
	    prip[2][nucpos] = 1;
	    secp[3][nucpos] = 1;
	    if (tr.traceACGT[0][tracePos] > tr.traceACGT[1][tracePos]) {
	      terp[0][nucpos] = 1;
	      quap[1][nucpos] = 1;
	    } else {
	      terp[1][nucpos] = 1;
	      quap[0][nucpos] = 1;
	    }
	  } else if ((pri[i] == 'T') && (sec[i] =='A')) {
	    prip[3][nucpos] = 1;
	    secp[0][nucpos] = 1;
	    if (tr.traceACGT[1][tracePos] > tr.traceACGT[2][tracePos]) {
	      terp[1][nucpos] = 1;
	      quap[2][nucpos] = 1;
	    } else {
	      terp[2][nucpos] = 1;
	      quap[1][nucpos] = 1;
	    }
	  } else if ((pri[i] == 'T') && (sec[i] =='C')) {
	    prip[3][nucpos] = 1;
	    secp[1][nucpos] = 1;
	    if (tr.traceACGT[0][tracePos] > tr.traceACGT[2][tracePos]) {
	      terp[0][nucpos] = 1;
	      quap[2][nucpos] = 1;
	    } else {
	      terp[2][nucpos] = 1;
	      quap[0][nucpos] = 1;
	    }
	  } else if ((pri[i] == 'T') && (sec[i] =='G')) {
	    prip[3][nucpos] = 1;
	    secp[2][nucpos] = 1;
	    if (tr.traceACGT[0][tracePos] > tr.traceACGT[1][tracePos]) {
	      terp[0][nucpos] = 1;
	      quap[1][nucpos] = 1;
	    } else {
	      terp[1][nucpos] = 1;
	      quap[0][nucpos] = 1;
	    }
	  }
	  ++nucpos;
	}	  
      }

      // Debug
      //for(uint32_t m = 0; m<4; ++m) {
      //for(uint32_t n = 0; n<diffnuc; ++n) {
      //  std::cout << tp[m][n] << ',';
      //}
      //std::cout << std::endl;
      //}
      
      // This is all tiny just do it brute-force
      for(uint32_t m = 0; m<4; ++m) {
	for(uint32_t n = 0; n<diffnuc; ++n) {
	  double pred = bestI * prip[m][n] + bestJ * secp[m][n] + bestK * terp[m][n] + bestL * quap[m][n];
	  bestSSE += (pred - tp[m][n]) * (pred - tp[m][n]);
	}
      }
      for(double i = 0; i <= 1; i += 0.01) {
	for(double j = 0; j <= 1; j += 0.01) {
	  if (i + j <= 1) {
	    for(double k = 0; k <= 1; k += 0.01) {
	      if (i + j + k <= 1) {
		double l = 1 - (i + j + k);
		double sse = 0;
		for(uint32_t m = 0; m<4; ++m) {
		  for(uint32_t n = 0; n<diffnuc; ++n) {
		    double pred = i * prip[m][n] + j * secp[m][n] + k * terp[m][n] + l * quap[m][n];
		    sse += (pred - tp[m][n]) * (pred - tp[m][n]);
		    if (sse >= bestSSE) break;
		  }
		}
		if (sse < bestSSE) {
		  bestSSE = sse;
		  bestL = l;
		  bestK = k;
		  bestJ = j;
		  bestI = i;
		}
	      }
	    }
	  }
	}
      }
    }
    //std::cout << bestSSE << ',' << bestI << ',' << bestJ << ',' << bestK << ',' << bestL << std::endl;
    return std::make_pair(bestI, bestJ);
  }
    
  

  template<typename TConfig, typename TDecomp>
  inline void
  writeDecomposition(TConfig const& c, TDecomp const& dcp) {
    boost::filesystem::path outdecomp(c.outfile.string() + ".decomp");
    std::ofstream ofile(outdecomp.string().c_str());
    ofile << "indel\tdecomp" << std::endl;
    for(uint32_t i = 0; i < dcp.size(); ++i) ofile << dcp[i].first << "\t" << dcp[i].second << std::endl;
    ofile.close();
  }  
}



#endif
