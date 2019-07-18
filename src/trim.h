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

#ifndef TRIM_H
#define TRIM_H

#include "abif.h"
#include "scf.h"
#include "json.h"
#include "fasta.h"

namespace tracy {

  template<typename TConfig>
  inline uint32_t
  nearestSNP(TConfig const& c, BaseCalls const& bc, uint32_t const rtp) {
    bool deadEnd = false;
    uint32_t offset = 0;
    while (!deadEnd) {
      deadEnd = true;
      if ((rtp + offset + c.trimRight < bc.secondary.size()) && (rtp + offset + c.trimRight < bc.primary.size())) {
	if (c.trimLeft < rtp + offset) {
	  if (bc.primary[rtp+offset] != bc.secondary[rtp+offset]) return (rtp + offset - c.trimLeft);
	}
	deadEnd = false;
      }
      if (offset + c.trimLeft < rtp) {
	if (bc.primary[rtp-offset] != bc.secondary[rtp-offset]) return (rtp - offset - c.trimLeft);
	deadEnd = false;
      }
      ++offset;
    }
    // No SNP found
    if (rtp > c.trimLeft) return rtp - c.trimLeft;
    else return c.trimLeft; 
  }
  
  template<typename TConfig>
  inline void
  trimTrace(TConfig const& c, BaseCalls const& bc, uint32_t& leftTrim, uint32_t& rightTrim) {
    // Screening window
    uint32_t win = 10;
    std::vector<int32_t> penalty(bc.secondary.size(), 0);
    typedef std::pair<uint32_t, double> TIdxVal;
    TIdxVal idxval = findBestTraceSection(bc, penalty, win);
    uint32_t bestIdx = idxval.first;
    double perBasePenalty = c.trimStringency * idxval.second;

    // Walk outwards to estimate trim position
    rightTrim = bc.secondary.size();
    leftTrim = 0;
    double localPenalty = 0;
    for(uint32_t i = bestIdx; ((i < bestIdx + win) && (i < bc.secondary.size())); ++i) localPenalty += penalty[i];
    for(uint32_t i = bestIdx; ((i + win) < bc.secondary.size()); ++i) {
      localPenalty -= penalty[i];
      localPenalty += penalty[i+win];
      if (localPenalty > (perBasePenalty * win)) {
	rightTrim = i;
	break;
      }
    }
    localPenalty = 0;
    for(uint32_t i = bestIdx; ((i < bestIdx + win) && (i < bc.secondary.size())); ++i) localPenalty += penalty[i];
    int32_t i = bestIdx - 1;
    while (i >= 0) {
      localPenalty -= penalty[i + win];
      localPenalty += penalty[i];
      if (localPenalty > (perBasePenalty * win)) {
	leftTrim = i + win - 1;
	break;
      }
      --i;
    }
    if (rightTrim < bc.secondary.size()) rightTrim = bc.secondary.size() - rightTrim;
    else rightTrim = 0;
  }


  // Trim basecalls but keep the entire trace (all sampling positions)
  inline void
  trimTrace(Trace const& tr, BaseCalls const& bc, uint32_t const trimLeft, uint32_t const trimRight, BaseCalls& nbc) {
    typedef Trace::TValue TValue;
    
    // Last basecall
    uint32_t len = bc.primary.size() - trimRight;
    
    // Rewrite arrays
    uint32_t bcpos = 0;
    TValue idx = bc.bcPos[0];
    for(TValue tracePos = 0; tracePos < (TValue) tr.traceACGT[0].size(); ++tracePos) {
      if (idx == tracePos) {
	if ((bcpos >= trimLeft) && (bcpos < len)) {
	  nbc.bcPos.push_back(tracePos);
	  nbc.primary.push_back(bc.primary[bcpos]);
	  nbc.secondary.push_back(bc.secondary[bcpos]);
	  nbc.consensus.push_back(bc.consensus[bcpos]);
	}
	if (bcpos < bc.bcPos.size() - 1) idx = bc.bcPos[++bcpos];
      }
    }
  }


  inline char
  reverseComplement(char const c) {
    switch (c) {
    case 'A': return 'T';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'T': return 'A';
    case 'N': return 'N';
    case 'H': return 'D';
    case 'V': return 'B';
    case 'M': return 'K';
    case 'Y': return 'R';
    case 'D': return 'H';
    case 'B': return 'V';
    case 'K': return 'M';
    case 'R': return 'Y';
    case 'U': return 'A';
    case 'S': return 'S';
    case 'W': return 'W';
    default: return c;
    }
  }
  
  inline void
  reverseComplementTrace(Trace const& tr, BaseCalls const& bc, Trace& ntr, BaseCalls& nbc) {
    typedef Trace::TMountains TMountains; 
    typedef Trace::TValue TValue;

    // Rewrite arrays
    uint32_t bcpos = bc.bcPos.size()-1;
    TValue idx = bc.bcPos[bcpos];
    ntr.traceACGT.resize(4, TMountains());
    TValue newTracePos = 0;
    for(TValue tracePos = tr.traceACGT[0].size(); tracePos > 0; --tracePos, ++newTracePos) {
      if (idx == tracePos - 1) {
	nbc.bcPos.push_back(newTracePos);
	nbc.primary.push_back(reverseComplement(bc.primary[bcpos]));
	nbc.secondary.push_back(reverseComplement(bc.secondary[bcpos]));
	nbc.consensus.push_back(reverseComplement(bc.consensus[bcpos]));
	ntr.qual.push_back(tr.qual[bcpos]);
	if (bcpos > 0) idx = bc.bcPos[--bcpos];
      }
      ntr.traceACGT[3].push_back(tr.traceACGT[0][tracePos-1]);        
      ntr.traceACGT[2].push_back(tr.traceACGT[1][tracePos-1]);
      ntr.traceACGT[1].push_back(tr.traceACGT[2][tracePos-1]);
      ntr.traceACGT[0].push_back(tr.traceACGT[3][tracePos-1]);
    }
  }

}

#endif
