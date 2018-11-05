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
  inline void
  trimTrace(TConfig const& c, BaseCalls const& bc, uint32_t& leftTrim, uint32_t& rightTrim) {
    std::vector<int32_t> penalty(bc.secondary.size(), 0);

    // Screening window
    uint32_t win = 10;
    
    // Secondary basecalls != [ACGT]
    uint32_t halfwin = (uint32_t) (win / 2);
    int32_t ambiguous = 0;
    for (uint32_t i = 0; ((i < win) && (i < bc.secondary.size())); ++i) {
      if (isAmbiguous(bc.secondary[i])) ++ambiguous;
    }
    for(uint32_t i = 0; ((i < halfwin) && (i < bc.secondary.size())); ++i) penalty[i] = ambiguous;
    for(uint32_t i = win; i < bc.secondary.size(); ++i) {
      if (isAmbiguous(bc.secondary[i-win])) --ambiguous;
      if (isAmbiguous(bc.secondary[i])) ++ambiguous;
      penalty[i - halfwin] = ambiguous;
    }
    for(uint32_t i = bc.secondary.size() - halfwin; i < bc.secondary.size(); ++i) penalty[i] = ambiguous;

    // Mean basecall distance
    double meanDist = 0;
    for(uint32_t i = 1; i < bc.secondary.size(); ++i) meanDist += (bc.bcPos[i] - bc.bcPos[i-1]);
    meanDist /= (bc.secondary.size() - 1);
    
    // Peak distance
    uint32_t peakVar = 0;
    for(uint32_t i = 0; (i + win < bc.secondary.size()); ++i) {
      uint32_t oldPos = 0;
      if (i>0) oldPos = bc.bcPos[i-1];
      uint32_t minDist = bc.bcPos[bc.secondary.size() - 1];
      uint32_t maxDist = 0;
      for(uint32_t k = 0; k < win; ++k) {
	uint32_t dist = bc.bcPos[i+k] - oldPos;
	oldPos = bc.bcPos[i+k];
	if (dist < minDist) minDist = dist;
	if (dist > maxDist) maxDist = dist;
      }
      peakVar = (int32_t) ( (std::abs((double) maxDist - meanDist) + std::abs((double) minDist - meanDist)) / 2);
      penalty[i + halfwin] += peakVar;
      if (i == 0) {
	for(uint32_t k = 0; k < halfwin; ++k) penalty[k] += peakVar;
      }
    }
    for(uint32_t i = bc.secondary.size() - halfwin; i < bc.secondary.size(); ++i) penalty[i] += peakVar;

    // Try to identify best 10% window
    uint32_t sourcewin = (int32_t) (0.1 * bc.secondary.size());
    uint32_t bestIdx = 0;
    int32_t bestVal = 99999999;
    for(uint32_t i = 0; ((i + sourcewin) < bc.secondary.size()); ++i) {
      int32_t penval = 0;
      for(uint32_t k = 0; k < sourcewin; ++k) penval += penalty[i+k];
      if (penval < bestVal) {
	bestVal = penval;
	bestIdx = i + (int32_t) (sourcewin / 2);
      }
    }

    // Walk outwards to estimate trim position
    rightTrim = bc.secondary.size();
    leftTrim = 0;
    double perBasePenalty = c.trimStringency * ((double) bestVal / (double) sourcewin);
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
