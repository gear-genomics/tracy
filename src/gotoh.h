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

#ifndef GOTOH_H
#define GOTOH_H

#include <iostream>
#include "align.h"

namespace tracy
{

  template<typename TAlign, typename TAlignConfig, typename TScoreObject>
  inline int
  gotohString(std::string const& s1, std::string const& s2, TAlign& align, TAlignConfig const& ac, TScoreObject const& sc)
  {
    typedef typename TScoreObject::TValue TScoreValue;

    // DP Matrix
    typedef boost::multi_array<TScoreValue, 2> TMatrix;
    std::size_t m = s1.size();
    std::size_t n = s2.size();
    TMatrix s(boost::extents[m+1][n+1]);
    TMatrix h(boost::extents[m+1][n+1]);
    TMatrix v(boost::extents[m+1][n+1]);

    // Initialization
    for(std::size_t col = 1; col <= n; ++col) {
      v[0][col] = -sc.inf;
      s[0][col] = _horizontalGap(ac, 0, m, sc.go + col * sc.ge);
      h[0][col] = _horizontalGap(ac, 0, m, sc.go + col * sc.ge);
    }
    for(std::size_t row = 1; row <= m; ++row) {
      h[row][0] = -sc.inf;
      s[row][0] = _verticalGap(ac, 0, n, sc.go + row * sc.ge);
      v[row][0] = _verticalGap(ac, 0, n, sc.go + row * sc.ge);
    }
    s[0][0] = 0;
    v[0][0] = -sc.inf;
    h[0][0] = -sc.inf;

    // Recursion
    for(std::size_t row = 1; row <= m; ++row) {
      for(std::size_t col = 1; col <= n; ++col) {
	h[row][col] = std::max(s[row][col-1] + _horizontalGap(ac, row, m, sc.go + sc.ge), h[row][col-1] + _horizontalGap(ac, row, m, sc.ge));
	v[row][col] = std::max(s[row-1][col] + _verticalGap(ac, col, n, sc.go + sc.ge), v[row-1][col] + _verticalGap(ac, col, n, sc.ge));
	s[row][col] = std::max(std::max(s[row-1][col-1] + _scoreString(s1, s2, row-1, col-1, sc), h[row][col]), v[row][col]);
      }
    }

    // Trace-back
    std::size_t row = m;
    std::size_t col = n;
    char lastMatrix = 's';
    typedef std::vector<char> TTrace;
    TTrace trace;
    while ((row>0) || (col>0)) {
      if (lastMatrix == 's') {
	if (s[row][col] == h[row][col]) lastMatrix = 'h';
	else if (s[row][col] == v[row][col]) lastMatrix = 'v';
	else {
	  --row;
	  --col;
	  trace.push_back('s');
	}
      } else if (lastMatrix == 'h') {
	if (h[row][col] != h[row][col-1] + _horizontalGap(ac, row, m, sc.ge)) lastMatrix = 's';
	--col;
	trace.push_back('h');
      } else if (lastMatrix == 'v') {
	if (v[row][col] != v[row-1][col] + _verticalGap(ac, col, n, sc.ge)) lastMatrix = 's';
	--row;
	trace.push_back('v');
      }
    }

    // Create alignment
    _createAlignmentString(trace, s1, s2, align);

    // Score
    return s[m][n];
  }


  template<typename TChar, typename TProfile>
  inline void
  _createMSAProfile(boost::multi_array<TChar, 2> const& a, TProfile& p)
  {
    typedef typename boost::multi_array<TChar, 2>::index TAIndex;
    typedef typename TProfile::index TPIndex;
    p.resize(boost::extents[6][a.shape()[1]]);   // 'A', 'C', 'G', 'T', 'N', '-'
    for (TAIndex j = 0; j < (TAIndex) a.shape()[1]; ++j) {
      for(TPIndex k = 0; k < 6; ++k) p[k][j] = 0;
      int sum = 0;
      for(TAIndex i = 0; i < (TAIndex) a.shape()[0]; ++i) {
	++sum;
	if ((a[i][j] == 'A') || (a[i][j] == 'a')) p[0][j] += 1;
	else if ((a[i][j] == 'C') || (a[i][j] == 'c')) p[1][j] += 1;
	else if ((a[i][j] == 'G') || (a[i][j] == 'g')) p[2][j] += 1;
	else if ((a[i][j] == 'T') || (a[i][j] == 't')) p[3][j] += 1;
	else if ((a[i][j] == 'N') || (a[i][j] == 'n')) p[4][j] += 1;
	else if (a[i][j] == '-') p[5][j] += 1;
	else --sum;
      }
      for(TPIndex k = 0; k<6; ++k) p[k][j] /= sum;
    }
  }

  template<typename TTrace, typename TChar, typename TAlign>
  inline void
  _createMSAAlignment(TTrace const& trace, boost::multi_array<TChar, 2> const& a1, boost::multi_array<TChar, 2> const& a2, TAlign& align)
  {
    typedef typename TAlign::index TAIndex;
    TAIndex numN = a1.shape()[0];
    TAIndex numM = a2.shape()[0];
    align.resize(boost::extents[numN + numM][trace.size()]);
    TAIndex row = 0;
    TAIndex col = 0;
    TAIndex ai = 0;
    for(typename TTrace::const_reverse_iterator itT = trace.rbegin(); itT != trace.rend(); ++itT, ++ai) {
      if (*itT == 's') {
	for(TAIndex i = 0; i<numN; ++i) align[i][ai] = a1[i][row];
	for(TAIndex i = 0; i<numM; ++i) align[numN + i][ai] = a2[i][col];
	++row;
	++col;
      } else if (*itT =='h') {
	for(TAIndex i = 0; i<numN; ++i) align[i][ai] = '-';
	for(TAIndex i = 0; i<numM; ++i) align[numN + i][ai] = a2[i][col];
	++col;
      } else {
	for(TAIndex i = 0; i<numN; ++i) align[i][ai] = a1[i][row];
	for(TAIndex i = 0; i<numM; ++i) align[numN + i][ai] = '-';
	++row;
      }
    }
  }
  
  template<typename TAlign1, typename TAlign2, typename TAlign, typename TAlignConfig, typename TScoreObject>
  inline int
  gotohMSA(TAlign1 const& a1, TAlign2 const& a2, TAlign& align, TAlignConfig const& ac, TScoreObject const& sc)
  {
    typedef typename TScoreObject::TValue TScoreValue;

    // DP Matrix
    typedef boost::multi_array<TScoreValue, 2> TMatrix;
    std::size_t m = a1.shape()[1];
    std::size_t n = a2.shape()[1];
    TMatrix s(boost::extents[m+1][n+1]);
    TMatrix h(boost::extents[m+1][n+1]);
    TMatrix v(boost::extents[m+1][n+1]);

    // Create profile
    typedef boost::multi_array<double, 2> TProfile;
    TProfile p1;
    _createMSAProfile(a1, p1);
    TProfile p2;
    _createMSAProfile(a2, p2);

    // Initialization
    for(std::size_t col = 1; col <= n; ++col) {
      v[0][col] = -sc.inf;
      s[0][col] = _horizontalGap(ac, 0, m, sc.go + col * sc.ge);
      h[0][col] = _horizontalGap(ac, 0, m, sc.go + col * sc.ge);
    }
    for(std::size_t row = 1; row <= m; ++row) {
      h[row][0] = -sc.inf;
      s[row][0] = _verticalGap(ac, 0, n, sc.go + row * sc.ge);
      v[row][0] = _verticalGap(ac, 0, n, sc.go + row * sc.ge);
    }
    s[0][0] = 0;
    v[0][0] = -sc.inf;
    h[0][0] = -sc.inf;

    // Recursion
    for(std::size_t row = 1; row <= m; ++row) {
      for(std::size_t col = 1; col <= n; ++col) {
	h[row][col] = std::max(s[row][col-1] + _horizontalGap(ac, row, m, sc.go + sc.ge), h[row][col-1] + _horizontalGap(ac, row, m, sc.ge));
	v[row][col] = std::max(s[row-1][col] + _verticalGap(ac, col, n, sc.go + sc.ge), v[row-1][col] + _verticalGap(ac, col, n, sc.ge));
	s[row][col] = std::max(std::max(s[row-1][col-1] + _score(p1, p2, row-1, col-1, sc), h[row][col]), v[row][col]);
      }
    }

    // Trace-back
    std::size_t row = m;
    std::size_t col = n;
    char lastMatrix = 's';
    typedef std::vector<char> TTrace;
    TTrace trace;
    while ((row>0) || (col>0)) {
      if (lastMatrix == 's') {
	if (s[row][col] == h[row][col]) lastMatrix = 'h';
	else if (s[row][col] == v[row][col]) lastMatrix = 'v';
	else {
	  --row;
	  --col;
	  trace.push_back('s');
	  //std::cerr << a1[0][row] << a2[0][col] << std::endl;
	}
      } else if (lastMatrix == 'h') {
	if (h[row][col] != h[row][col-1] + _horizontalGap(ac, row, m, sc.ge)) lastMatrix = 's';
	--col;
	trace.push_back('h');
	//std::cerr << '-' << a2[0][col] << std::endl;
      } else if (lastMatrix == 'v') {
	if (v[row][col] != v[row-1][col] + _verticalGap(ac, col, n, sc.ge)) lastMatrix = 's';
	--row;
	trace.push_back('v');
	//std::cerr << a1[0][row] << '-' << std::endl;
      }
    }

    // Create alignment
    _createMSAAlignment(trace, a1, a2, align);

    // Score
    return s[m][n];
  }

  
  
  template<typename TProfile, typename TAlign, typename TAlignConfig, typename TScoreObject>
  inline int
  gotoh(TProfile const& p1, TProfile const& p2, TAlign& align, TAlignConfig const& ac, TScoreObject const& sc) {
    typedef typename TScoreObject::TValue TScoreValue;

    // DP Matrix
    typedef boost::multi_array<TScoreValue, 2> TMatrix;
    std::size_t m = p1.shape()[1];
    std::size_t n = p2.shape()[1];
    TMatrix s(boost::extents[m+1][n+1]);
    TMatrix h(boost::extents[m+1][n+1]);
    TMatrix v(boost::extents[m+1][n+1]);

    // Initialization
    for(std::size_t col = 1; col <= n; ++col) {
      v[0][col] = -sc.inf;
      s[0][col] = _horizontalGap(ac, 0, m, sc.go + col * sc.ge);
      h[0][col] = _horizontalGap(ac, 0, m, sc.go + col * sc.ge);
    }
    for(std::size_t row = 1; row <= m; ++row) {
      h[row][0] = -sc.inf;
      s[row][0] = _verticalGap(ac, 0, n, sc.go + row * sc.ge);
      v[row][0] = _verticalGap(ac, 0, n, sc.go + row * sc.ge);
    }
    s[0][0] = 0;
    v[0][0] = -sc.inf;
    h[0][0] = -sc.inf;

    // Recursion
    for(std::size_t row = 1; row <= m; ++row) {
      for(std::size_t col = 1; col <= n; ++col) {
	h[row][col] = std::max(s[row][col-1] + _horizontalGap(ac, row, m, sc.go + sc.ge), h[row][col-1] + _horizontalGap(ac, row, m, sc.ge));
	v[row][col] = std::max(s[row-1][col] + _verticalGap(ac, col, n, sc.go + sc.ge), v[row-1][col] + _verticalGap(ac, col, n, sc.ge));
	s[row][col] = std::max(std::max(s[row-1][col-1] + _score(p1, p2, row-1, col-1, sc), h[row][col]), v[row][col]);
      }
    }

    // Trace-back
    std::size_t row = m;
    std::size_t col = n;
    char lastMatrix = 's';
    typedef std::vector<char> TTrace;
    TTrace trace;
    while ((row>0) || (col>0)) {
      if (lastMatrix == 's') {
	if (s[row][col] == h[row][col]) lastMatrix = 'h';
	else if (s[row][col] == v[row][col]) lastMatrix = 'v';
	else {
	  --row;
	  --col;
	  trace.push_back('s');
	}
      } else if (lastMatrix == 'h') {
	if (h[row][col] != h[row][col-1] + _horizontalGap(ac, row, m, sc.ge)) lastMatrix = 's';
	--col;
	trace.push_back('h');
      } else if (lastMatrix == 'v') {
	if (v[row][col] != v[row-1][col] + _verticalGap(ac, col, n, sc.ge)) lastMatrix = 's';
	--row;
	trace.push_back('v');
      }
    }

    // Create alignment
    _createAlignment(trace, p1, p2, align);

    // Score
    return s[m][n];
  }

}

#endif
