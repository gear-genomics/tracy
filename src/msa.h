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

#ifndef MSA_H
#define MSA_H

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/multi_array.hpp>
#include "gotoh.h"
#include "needle.h"

namespace tracy {

  
  template<typename TArray>
  inline void
  overwriteArray(TArray const& in, TArray& out) {
    typedef typename TArray::index TAIndex;
    out.resize(boost::extents[in.shape()[0]][in.shape()[1]]);
    //std::cerr << "Sequence alignment" << std::endl;
    for(TAIndex i = 0; i < (TAIndex) in.shape()[0]; ++i) {
      for(TAIndex j = 0; j < (TAIndex) in.shape()[1]; ++j) {
	out[i][j] = in[i][j];
	//std::cerr << out[i][j];
      }
      //std::cerr << std::endl;
    }
    //std::cerr << std::endl;
  }

  
  template<typename TConfig, typename TSeqProfiles, typename TDistArray>
  inline void
  distanceMatrix(TConfig const& c, TSeqProfiles const& sps, TDistArray& d) {
    for (uint32_t i = 0; i < sps.size(); ++i) {
      for (uint32_t j = i+1; j < sps.size(); ++j) {
	AlignConfig<true, true> alignconf;
	d[i][j] = gotohScore(sps[i], sps[j], alignconf, c.aliscore);
      }
    }
  }

  template<typename TDistArray, typename TDIndex>
  inline int
  closestPair(TDistArray const& d, TDIndex num, TDIndex& dI, TDIndex& dJ) {
    int dMax = -1;
    for (TDIndex i = 0; i<num; ++i) {
      for (TDIndex j = i+1; j<num; ++j) {
	if (d[i][j]>dMax) {
	  dMax = d[i][j];
	  dI = i;
	  dJ = j;
	}
      }
    }
    return dMax;
  }

  template<typename TDistArray, typename TPhylogeny, typename TDIndex>
  inline void
  updateDistanceMatrix(TDistArray& d, TPhylogeny const& p, TDIndex num, TDIndex& dI, TDIndex& dJ) {
    for (TDIndex i = 0; i < num; ++i) 
      if (p[i][0] == -1) 
	d[i][num] = (((dI < i) ? d[dI][i] : d[i][dI]) + ((dJ < i) ? d[dJ][i] : d[i][dJ])) / 2;
    for (TDIndex i = 0; i<dI; ++i) d[i][dI] = -1;
    for (TDIndex i = dI+1; i<num+1; ++i) d[dI][i] = -1;
    for (TDIndex i = 0; i<dJ; ++i) d[i][dJ] = -1;
    for (TDIndex i = dJ+1; i<num+1; ++i) d[dJ][i] = -1;
  }

  template<typename TDistArray, typename TPhylogeny, typename TDIndex>
  inline TDIndex
  upgma(TDistArray& d, TPhylogeny& p, TDIndex num) {
    TDIndex nn = num;
    for(;nn<2*num+1; ++nn) {
      TDIndex dI = 0;
      TDIndex dJ = 0;
      if (closestPair(d, nn, dI, dJ) == -1) break;
      p[dI][0] = nn;
      p[dJ][0] = nn;
      p[nn][1] = dI;
      p[nn][2] = dJ;
      updateDistanceMatrix(d, p, nn, dI, dJ);
    }
    return (nn > 0) ? (nn - 1) : 0;
  }

  template<typename TConfig, typename TSeqProfiles, typename TPhylogeny, typename TDIndex, typename TAlign, typename TProfile, typename TSeqIdx>
  inline void
  palign(TConfig const& c, TSeqProfiles const& sps, TPhylogeny const& p, TDIndex root, TAlign& align, TProfile& prof, TSeqIdx& sidx) {
    if ((p[root][1] == -1) && (p[root][2] == -1)) {
      align.resize(boost::extents[1][sps[root].shape()[1]]);
      for(uint32_t ind = 0; ind < sps[root].shape()[1]; ++ind) align[0][ind] = _profileConsChar(sps[root], ind);
      copyProfile(sps[root], prof);
      sidx.push_back(root);
    } else {
      TAlign align1;
      TProfile prof1;
      TSeqIdx sidx1;
      palign(c, sps, p, p[root][1], align1, prof1, sidx1);
      TAlign align2;
      TProfile prof2;
      TSeqIdx sidx2;
      palign(c, sps, p, p[root][2], align2, prof2, sidx2);
      AlignConfig<true, true> endFreeAlign;

      // Debug
      //std::cerr << prof1.shape()[0] << ',' << prof1.shape()[1] << std::endl;
      //std::cerr << align1.shape()[0] << ',' << align1.shape()[1] << std::endl;
      //std::cerr << prof2.shape()[0] << ',' << prof2.shape()[1] << std::endl;
      //std::cerr << align2.shape()[0] << ',' << align2.shape()[1] << std::endl;
      
      // Profile-to-profile alignment
      TAlign alignNew;
      gotoh(prof1, prof2, alignNew, endFreeAlign, c.aliscore);

      // Debug profile alignment
      //std::cerr << "Profile alignment" << std::endl;
      //for(uint32_t i = 0; i < alignNew.shape()[0]; ++i) {
      //for(uint32_t j = 0; j < alignNew.shape()[1]; ++j) std::cerr << alignNew[i][j];
	//std::cerr << std::endl;
      //}
      
      // Create new sequence alignment based on profile alignment
      TAlign alignCombined;
      uint32_t nSeq = align1.shape()[0] + align2.shape()[0];
      uint32_t nCol = alignNew.shape()[1];
      uint32_t a1p = 0;
      uint32_t a2p = 0;
      alignCombined.resize(boost::extents[nSeq][nCol]);
      for(uint32_t j = 0; j < nCol; ++j) {
	if (alignNew[0][j] != '-') {
	  for(uint32_t k = 0; k < align1.shape()[0]; ++k) alignCombined[k][j] = align1[k][a1p];
	  ++a1p;
	} else {
	  for(uint32_t k = 0; k < align1.shape()[0]; ++k) alignCombined[k][j] = '-';
	}
	if (alignNew[1][j] != '-') {
	  uint32_t ind = 0;
	  for(uint32_t k = align1.shape()[0]; k < nSeq; ++k) alignCombined[k][j] = align2[ind++][a2p];
	  ++a2p;
	} else {
	  for(uint32_t k = align1.shape()[0]; k < nSeq; ++k) alignCombined[k][j] = '-';
	}
      }

      // Overwrite old alignment                                                   
      overwriteArray(alignCombined, align);

      // Create alignment profile
      _createProfile(align, prof);

      // Create sequence index
      sidx.resize(nSeq);
      for(uint32_t k = 0; k < sidx1.size(); ++k) sidx[k] = sidx1[k];
      uint32_t ind = 0;
      for(uint32_t k = sidx1.size(); k < nSeq; ++k) sidx[k] = sidx2[ind++];
    }
  }

  template<typename TConfig, typename TAlign>
  inline void
  consensus(TConfig const& c, TAlign const& align, std::string& gapped, std::string& cs, bool ignoreLast) {
    typedef typename TAlign::index TAIndex;

    // Ignore last sequence?
    uint32_t getRidOffRef = 0;
    if (ignoreLast) getRidOffRef = 1;
    
    // Calculate coverage
    typedef boost::multi_array<bool, 2> TFlag;
    TFlag fl;
    fl.resize(boost::extents[align.shape()[0] - getRidOffRef][align.shape()[1]]);
    typedef std::vector<int32_t> TCoverage;
    TCoverage cov;
    cov.resize(align.shape()[1], 0);
    for(TAIndex i = 0; (i < ((TAIndex) align.shape()[0] - getRidOffRef)); ++i) {
      int start = 0;
      int end = -1;
      for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
	fl[i][j] = false;
	if (align[i][j] != '-') end = j;
	else if (end == -1) start = j + 1;
      }
      for(TAIndex j = start; j<=end; ++j) {
	++cov[j];
	fl[i][j] = true;
      }
    }

    // Minimum number of aligned sequences
    int32_t covThreshold = (int32_t) (c.fractionCalled * (align.shape()[0] - getRidOffRef));
    TAIndex j = 0;
    std::vector<char> cons(align.shape()[1], '-');
    for(typename TCoverage::const_iterator itCov = cov.begin(); itCov != cov.end(); ++itCov, ++j) {
      int32_t maxIdx = 4;  // Leading/trailing gaps until min. coverage is reached
      if ((*itCov >= 1) && (*itCov >= covThreshold)) {
	// Get consensus letter
	std::vector<int32_t> count(5, 0); // ACGT-
	for(TAIndex i = 0; (i < ((TAIndex) align.shape()[0] - getRidOffRef)); ++i) {
	  if (fl[i][j]) {
	    if ((align[i][j] == 'A') || (align[i][j] == 'a')) ++count[0];
	    else if ((align[i][j] == 'C') || (align[i][j] == 'c')) ++count[1];
	    else if ((align[i][j] == 'G') || (align[i][j] == 'g')) ++count[2];
	    else if ((align[i][j] == 'T') || (align[i][j] == 't')) ++count[3];
	    else ++count[4];
	  }
	}
	maxIdx = 0;
	int32_t maxCount = count[0];
	for(uint32_t i = 1; i<5; ++i) {
	  if (count[i] > maxCount) {
	    maxCount = count[i];
	    maxIdx = i;
	  }
	}
      }
      switch (maxIdx) {
      case 0: cons[j] = 'A'; break;
      case 1: cons[j] = 'C'; break;
      case 2: cons[j] = 'G'; break;
      case 3: cons[j] = 'T'; break;
      default: break;
      }
    }
    gapped = std::string(cons.begin(), cons.end());
    for(uint32_t i = 0; i<cons.size(); ++i) {
      if (cons[i] != '-') cs.push_back(cons[i]);
    }
  }



  template<typename TConfig, typename TSeqProfiles>
  inline void
  revSeqBasedOnDist(TConfig const& c, TSeqProfiles& seq, std::vector<bool>& fwd) {
    typedef typename TSeqProfiles::value_type TProfile;
    
    // Compute distance matrix
    int32_t totalScore = 0;
    typedef boost::multi_array<int32_t, 2> TDistArray;
    typedef typename TDistArray::index TDIndex;
    TDIndex num = seq.size();
    TDistArray d(boost::extents[num][num]);
    for (TDIndex i = 0; i<num; ++i) {
      d[i][i] = 0;
      for (TDIndex j = i+1; j<num; ++j) {
	AlignConfig<true, true> alignconf;
	d[i][j] = gotohScore(seq[i], seq[j], alignconf, c.aliscore);
	d[j][i] = d[i][j];
	totalScore += d[i][j];
      }
    }
    // Debug: Initial score
    //std::cerr << "Score: " << totalScore << std::endl;

    // Optimize fwd-rev X-times
    bool iterateScore = true;
    int32_t updatedScore = 0;
    while (iterateScore) {
      // Get quality
      typedef std::pair<int32_t, int32_t> TSeqScore;
      std::vector<TSeqScore> seqQuality;
      for (TDIndex i = 0; i<num; ++i) {
	int32_t rowSum = 0;
	for (TDIndex j = 0; j<num; ++j) {
	  rowSum += d[i][j];
	}
	seqQuality.push_back(std::make_pair(rowSum, i));
      }
      // Sort by worst sequence
      std::sort(seqQuality.begin(), seqQuality.end());

      // Update scores
      for(uint32_t k = 0; k<seqQuality.size(); ++k) {
	TProfile s;
        reverseComplementProfile(seq[seqQuality[k].second], s);
	std::vector<int32_t> newD(num, 0);
	int32_t scoreSum = 0;
	int32_t oldScoreSum = 0;
	for (TDIndex i = 0; i<num; ++i) {
	  if (i != seqQuality[k].second) {
	    AlignConfig<true, true> alignconf;
	    newD[i] = gotohScore(seq[i], s, alignconf, c.aliscore);
	    oldScoreSum += d[i][seqQuality[k].second];
	    scoreSum += newD[i];
	  }
	}
	if (scoreSum >= oldScoreSum) {
	  seq[seqQuality[k].second] = s;
	  fwd[seqQuality[k].second] = (!fwd[seqQuality[k].second]);
	  for (TDIndex i = 0; i<num; ++i) {
	    d[i][seqQuality[k].second] = newD[i];
	    d[seqQuality[k].second][i] = d[i][seqQuality[k].second];
	  }
	}

	// Updated Score
	updatedScore = 0;
	for (TDIndex i = 0; i<num; ++i) {
	  for (TDIndex j = 0; j<num; ++j) {
	    updatedScore += d[i][j];
	  }
	}
	std::cout << "." << std::flush;
      }
      
      // Updated score
      updatedScore = 0;
      for (TDIndex i = 0; i<num; ++i) {
	for (TDIndex j = 0; j<num; ++j) {
	  updatedScore += d[i][j];
	}
      }
      if (totalScore < updatedScore) totalScore = updatedScore;
      else iterateScore = false;
    }
    std::cout << std::endl;
  }

  template<typename TConfig, typename TSeqProfiles, typename TAlign>
  inline void
  msa(TConfig const& c, TSeqProfiles const& sps, TAlign& align, std::vector<uint32_t>& seqidx) {
    typedef typename TSeqProfiles::value_type TProfile;
    
    // Compute distance matrix
    typedef boost::multi_array<int32_t, 2> TDistArray;
    typedef typename TDistArray::index TDIndex;
    TDIndex num = sps.size();
    TDistArray d(boost::extents[2*num+1][2*num+1]);
    for (TDIndex i = 0; i<(2*num+1); ++i) 
      for (TDIndex j = i+1; j<(2*num+1); ++j) 
	d[i][j]=-1;
    distanceMatrix(c, sps, d);

    // UPGMA
    typedef boost::multi_array<int, 2> TPhylogeny;
    TPhylogeny p(boost::extents[2*num+1][3]);
    for(TDIndex i = 0; i<(2*num+1); ++i) 
      for (TDIndex j = 0; j<3; ++j) p[i][j] = -1;
    TDIndex root = upgma(d, p, num);

    // Debug guide tree
    //std::cerr << "Phylogeny" << std::endl;
    //std::cerr << "#Sequences: " << sps.size() << std::endl;
    //std::cerr << "Root: " << root << std::endl;
    //std::cerr << "Node:Parent\tLeftChild\tRightChild" << std::endl;
    //for(TDIndex i = 0; i<(2*num+1); ++i) {
    //std::cerr << i << ':' << '\t';
    //for (TDIndex j = 0; j<3; ++j) {
    //std::cerr << p[i][j] << '\t';
    //}
    //std::cerr << std::endl;
    //}
    
    // Progressive Alignment
    TProfile prof;
    palign(c, sps, p, root, align, prof, seqidx);
  }

}

#endif
