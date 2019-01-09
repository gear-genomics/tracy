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

}

#endif
