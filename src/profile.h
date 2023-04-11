#ifndef PROFILE_H
#define PROFILE_H

namespace tracy
{

  inline bool
  _inBaseCalled(uint32_t const k, char const p, char const s) {
    if (k == 0) {
      if ((p == 'A') || (p == 'R') || (p == 'W') || (p == 'M') || (s == 'A') || (s == 'R') || (s == 'W') || (s == 'M')) return true;
    } else if (k == 1) {
      if ((p == 'C') || (p == 'Y') || (p == 'S') || (p == 'M') || (s == 'C') || (s == 'Y') || (s == 'S') || (s == 'M')) return true;
    } else if (k == 2) {
      if ((p == 'G') || (p == 'R') || (p == 'S') || (p == 'K') || (s == 'G') || (s == 'R') || (s == 'S') || (s == 'K')) return true;
    } else if (k == 3) {
      if ((p == 'T') || (p == 'Y') || (p == 'W') || (p == 'K') || (s == 'T') || (s == 'Y') || (s == 'W') || (s == 'K')) return true;
    }
    return false;
  }

  template<typename TProfile>
  inline void
  createProfile(Trace const& tr, BaseCalls const& bc, TProfile& p, int32_t trimleft, int32_t trimright) {
    if (trimleft + trimright >= (int32_t) bc.bcPos.size()) {
      trimleft = 0;
      trimright = 0;
    }
    int32_t sz = bc.bcPos.size() - (trimleft + trimright);
    p.resize(boost::extents[6][sz]);   // 'A', 'C', 'G', 'T', 'N', '-'
    for(int32_t j = trimleft; j < (trimleft + sz); ++j) {
      float totalsig = 0;
      float allBaseSig = 0;
      for(uint32_t k = 0; k<4; ++k) {
	allBaseSig += tr.traceACGT[k][bc.bcPos[j]];
	if (_inBaseCalled(k, bc.primary[j], bc.secondary[j])) totalsig += tr.traceACGT[k][bc.bcPos[j]];
      }
      p[4][j-trimleft] = 0;
      p[5][j-trimleft] = 0;
      if (totalsig == 0) {
	for(uint32_t k = 0; k<4; ++k) p[k][j-trimleft] = 0.25;
      } else {
	for(uint32_t k = 0; k<4; ++k) {
	  if (_inBaseCalled(k, bc.primary[j], bc.secondary[j])) p[k][j-trimleft] = ((float) (tr.traceACGT[k][bc.bcPos[j]]) / totalsig);
	}
	// Sometimes basecall signal is tiny fraction of all base signals (e.g., missing peaks in signal ramps)
	float normfac = totalsig / allBaseSig;
	for(uint32_t k = 0; k<4; ++k) {
	  p[k][j-trimleft] = normfac * p[k][j-trimleft] + (1 - normfac) * 0.25;
	}
      }
    }
  }

  template<typename TProfile>
  inline void
  createProfile(Trace const& tr, BaseCalls const& bc, TProfile& p) {
    createProfile(tr, bc, p, 0, 0);
  }
  
  template<typename TConfig, typename TProfile>
  inline void
  createProfile(TConfig const& c, ReferenceSlice const& rs, TProfile& p) {
    if (rs.filetype != 2) return _createProfile(rs.refslice, p);
    else {
      // Create profile from trace
      Trace wt;
      if (!readab(c.genome.string(), wt)) return;
      BaseCalls wtbc;
      basecall(wt, wtbc, c.pratio);
      return createProfile(wt, wtbc, p);
    }
  }

  template<typename TProfile>
  inline void
  reverseComplementProfile(TProfile const& p, TProfile& out) {
    out.resize(boost::extents[6][p.shape()[1]]);   // 'A', 'C', 'G', 'T', 'N', '-'
    int32_t pIdx = p.shape()[1] - 1;
    int32_t outIdx = 0;
    while (pIdx >= 0) {
      out[0][outIdx] = p[3][pIdx];
      out[1][outIdx] = p[2][pIdx];
      out[2][outIdx] = p[1][pIdx];
      out[3][outIdx] = p[0][pIdx];
      out[4][outIdx] = p[4][pIdx];
      out[5][outIdx] = p[5][pIdx];
      --pIdx;
      ++outIdx;
    }
  }

  template<typename TProfile>
  inline void
  copyProfile(TProfile const& p, TProfile& out) {
    out.resize(boost::extents[p.shape()[0]][p.shape()[1]]);
    for(uint32_t i = 0; i < p.shape()[0]; ++i) {
      for(uint32_t j = 0; j < p.shape()[1]; ++j) out[i][j] = p[i][j];
    }
  }
}

#endif
