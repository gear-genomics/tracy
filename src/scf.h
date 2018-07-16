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

#ifndef SCF_H
#define SCF_H

#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <string>
#include <iterator>
#include <algorithm>
#include <stdint.h>

#include "abif.h"

namespace tracy
{

inline int32_t
traceFormat(std::string const& filename) {
  // Read the mountains
  std::ifstream bfile(filename.c_str(), std::ios_base::binary | std::ios::ate);
  std::streamsize bsize = bfile.tellg();
  bfile.seekg(0, std::ios::beg);
  std::vector<char> buffer(bsize);
  int32_t ft = -1;
  if (bfile.read(buffer.data(), bsize)) {
    std::string filetype = readBinStr(buffer, 0, 4);
    if (filetype == "ABIF") ft = 0;
    else if (filetype == ".scf") ft = 1;
  }
  bfile.close();
  return ft;
}

  
  
inline bool
readscf(std::string const& filename, Trace& tr) {
  typedef Trace::TACGTMountains TACGTMountains;
  typedef TACGTMountains::value_type TMountains;
  tr.traceACGT.resize(4, TMountains());
    
  // Read the mountains
  std::ifstream bfile(filename.c_str(), std::ios_base::binary | std::ios::ate);
  std::streamsize bsize = bfile.tellg();
  bfile.seekg(0, std::ios::beg);
  std::vector<char> buffer(bsize);
  if (bfile.read(buffer.data(), bsize)) {
    std::string filetype = readBinStr(buffer, 0, 4);
    if (filetype != ".scf") {
      std::cerr << "File is not in SCF format!" << std::endl;
      bfile.close();
      return false;
    }

    // Header
    int32_t numSamplings = readBinI32(buffer, 4);
    int32_t offset = readBinI32(buffer, 8);
    int32_t numBases = readBinI32(buffer, 12);
    int32_t basesOffset = readBinI32(buffer, 24);
    std::string version = readBinStr(buffer, 36, 4);
    float numericVersion = boost::lexical_cast<float>(version);

    // Trace signal
    if (numericVersion > 2.9) {
      for(int32_t i = 0; i < 4; ++i) {
	for(int32_t k = i * numSamplings; k < (i + 1) * numSamplings; ++k) {
	  tr.traceACGT[i].push_back(readBinI16(buffer, offset + k * 2));
	}
	// offset
	for(uint32_t k = 0; k < 2; ++k) {
	  int16_t prev = 0;
	  for(uint32_t p = 0; p<tr.traceACGT[i].size(); ++p) {
	    tr.traceACGT[i][p] += prev;
	    prev = tr.traceACGT[i][p];
	  }
	}
	// Debug
	//for(uint32_t p = 0; p<tr.traceACGT[i].size(); ++p) std::cerr << p << ": " << tr.traceACGT[i][p] << std::endl;
      }
    } else {
      for(int32_t k = 0; k < 4 * numSamplings; ++k) tr.traceACGT[k%4].push_back(readBinI16(buffer, offset + k * 2));
    }

    // Basecall position
    if (numericVersion > 2.9) {
      for(int32_t k = 0; k < numBases; ++k) {
	tr.basecallpos.push_back(readBinI32(buffer, basesOffset + k * 4));
	tr.qual.push_back(0);
      }
    } else {
      std::cerr << "SCF version greater 2.9 required!" << std::endl;
      return false;
    }
	
  }
  bfile.close();

  // Close input file
  return true;
}


}

#endif
