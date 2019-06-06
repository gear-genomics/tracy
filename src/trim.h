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


  inline void
  trimTrace(Trace const& tr, BaseCalls const& bc, uint32_t const trimLeft, uint32_t const trimRight, Trace& ntr, BaseCalls& nbc) {
    typedef Trace::TMountains TMountains; 
    typedef Trace::TValue TValue;

    // Last basecall
    uint32_t len = bc.primary.size() - trimRight;
    
    // Rewrite arrays
    uint32_t bcpos = 0;
    TValue idx = bc.bcPos[0];
    ntr.traceACGT.resize(4, TMountains());
    TValue leftOffset = tr.traceACGT[0].size();
    if (trimLeft == 0) leftOffset = 0;
    TValue rightOffset = tr.traceACGT[0].size();
    TValue newTracePos = 0;
    for(TValue tracePos = 0; tracePos < (TValue) tr.traceACGT[0].size(); ++tracePos) {
      if (idx == tracePos) {
	if ((tracePos >= leftOffset) && (tracePos <= rightOffset)) {
	  nbc.bcPos.push_back(newTracePos);
	  nbc.primary.push_back(bc.primary[bcpos]);
	  nbc.secondary.push_back(bc.secondary[bcpos]);
	  nbc.consensus.push_back(bc.consensus[bcpos]);
	  ntr.qual.push_back(tr.qual[bcpos]);
	}
	if (bcpos < bc.bcPos.size() - 1) idx = bc.bcPos[++bcpos];
	if (bcpos == trimLeft) {
	  TValue prevVal = 0;
	  if (bcpos > 0) prevVal = bc.bcPos[bcpos - 1];
	  leftOffset = (TValue) ((prevVal + idx) / 2);
	}
	if (bcpos + 1 == len) {
	  TValue nextVal = tr.traceACGT[0].size();
	  if (bcpos + 1 < bc.primary.size()) nextVal = bc.bcPos[bcpos + 1];
	  rightOffset = (TValue) ((idx + nextVal) / 2);
	}
      }
      if ((tracePos >= leftOffset) && (tracePos <= rightOffset)) {
	ntr.traceACGT[0].push_back(tr.traceACGT[0][tracePos]);        
	ntr.traceACGT[1].push_back(tr.traceACGT[1][tracePos]);
	ntr.traceACGT[2].push_back(tr.traceACGT[2][tracePos]);
	ntr.traceACGT[3].push_back(tr.traceACGT[3][tracePos]);
	++newTracePos;
      }
    }
  }

}

#endif
