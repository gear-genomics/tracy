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

#ifndef ABIF_H
#define ABIF_H

#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <string>
#include <iterator>
#include <algorithm>
#include <stdint.h>

namespace tracy
{

struct Abif {
  std::string key;
  std::string name;
  int32_t number;
  int16_t etype;
  int16_t esize;
  int32_t nelements;
  int32_t dsize;
  int32_t doffset;
};

struct Trace {
  typedef int32_t TValue;
  typedef std::vector<TValue> TMountains;
  typedef std::vector<TMountains> TACGTMountains;
  typedef std::vector<uint8_t> TQual;
  
  std::string basecalls1;
  std::string basecalls2;
  TQual qual;
  TMountains basecallpos;
  TACGTMountains traceACGT;
};


struct BaseCalls {
  typedef Trace::TValue TValue;
  typedef std::vector<TValue> TPosition;
  
  std::string consensus;
  std::string primary;
  std::string secondary;
  TPosition bcPos;
};


template<typename TValue>
inline std::string
_toString(TValue const a) {
  std::stringstream ss;
  ss << a;
  return ss.str();
}  
 
template<typename TACGTMountains, typename TMountains>
inline void
peak(TACGTMountains const& trace, float const s, float const e, TMountains& pVal, TMountains& pIdx) {
  typedef typename TMountains::value_type TValue;
  for(uint32_t k = 0; k<4; ++k) {
    TValue bestIdx = (int32_t) (std::floor(s)) + 1;
    TValue bestVal = (TValue) 0;
    for(int32_t i = std::max(1, (int32_t) std::floor(s)); i < std::min((int32_t) (trace[k].size() - 1), (int32_t) std::floor(e)); ++i) {
      if (((trace[k][i-1] <= trace[k][i]) && (trace[k][i] > trace[k][i+1])) || ((trace[k][i-1] < trace[k][i]) && (trace[k][i] >= trace[k][i+1]))) {
	if (trace[k][i] > bestVal) {
	  bestIdx = (TValue) i;
	  bestVal = (TValue) trace[k][i];
	}
      }
    }
    pVal.push_back(bestVal);
    pIdx.push_back(bestIdx);
  }
}

template<typename TMountains>
inline char
iupac(TMountains const& p) {
  if (p.size() == 1) {
    if (p[0] == 0) return 'A';
    else if (p[0] == 1) return 'C';
    else if (p[0] == 2) return 'G';
    else if (p[0] == 3) return 'T';
  } else if (p.size() == 2) {
    if ((p[0] == 0) && (p[1] == 2)) return 'R';
    else if ((p[0] == 1) && (p[1] == 3)) return 'Y';
    else if ((p[0] == 1) && (p[1] == 2)) return 'S';
    else if ((p[0] == 0) && (p[1] == 3)) return 'W';
    else if ((p[0] == 2) && (p[1] == 3)) return 'K';
    else if ((p[0] == 0) && (p[1])) return 'M';
  }
  return 'N';
}

inline char
iupac(char const one, char const two) {
  typedef Trace::TMountains TMountains;
  typedef TMountains::value_type TValue;
  TMountains p(2, 0);
  if (one == 'A') p[0] = 0;
  else if (one == 'C') p[0] = 1;
  else if (one == 'G') p[0] = 2;
  else if (one == 'T') p[0] = 3;
  if (two == 'A') p[1] = 0;
  else if (two == 'C') p[1] = 1;
  else if (two == 'G') p[1] = 2;
  else if (two == 'T') p[1] = 3;
  if (p[1] < p[0]) {
    TValue tmp = p[0];
    p[0] = p[1];
    p[1] = tmp;
  }
  return iupac(p);
}
  

inline std::string
readBinStr(std::vector<char> const& buffer, int32_t pos, int32_t len) {
  return std::string(buffer.begin() + pos, buffer.begin() + pos + len);
}

inline uint8_t
readBinUI8(std::vector<char> const& buffer, int32_t pos) {
  return (uint8_t)(buffer[pos]);
}

inline int32_t
readBinI32(std::vector<char> const& buffer, int32_t pos) {
  return (((uint32_t) 0) | ((uint8_t)(buffer[pos])<<24) | ((uint8_t)(buffer[pos+1])<<16) | ((uint8_t)(buffer[pos+2])<<8) | ((uint8_t)(buffer[pos+3])));
}

inline int16_t
readBinI16(std::vector<char> const& buffer, int32_t pos) {
  return (((uint16_t) 0) | ((uint8_t)(buffer[pos])<<8) | ((uint8_t)(buffer[pos+1])));
}

inline std::string
replaceNonDna(std::string const& str) {
  std::string out;
  for(uint32_t i = 0; i<str.size();++i) {
    if ((str[i] == 'A') || (str[i] == 'C') || (str[i] == 'G') || (str[i] == 'T')) out = out.append(str, i, 1);
    else out = out.append("N");
  }
  return out;
}

inline bool
readab(std::string const& filename, Trace& tr) {
  typedef Trace::TACGTMountains TACGTMountains;
  typedef TACGTMountains::value_type TMountains;
  TACGTMountains trace(4, TMountains());
  std::string acgtOrder;

  // Read the mountains
  std::ifstream bfile(filename.c_str(), std::ios_base::binary | std::ios::ate);
  std::streamsize bsize = bfile.tellg();
  bfile.seekg(0, std::ios::beg);
  std::vector<char> buffer(bsize);
  if (bfile.read(buffer.data(), bsize)) {
    std::string filetype = readBinStr(buffer, 0, 4);
    if (filetype != "ABIF") {
      std::cerr << "File is not in ABIF format!" << std::endl;
      bfile.close();
      return false;
    }
    //int16_t version = readBinI16(buffer, 4);
    //std::string name = readBinStr(buffer, 6, 4);
    //int32_t number = readBinI32(buffer, 10);
    //int16_t etype = readBinI16(buffer, 14);
    int16_t esize = readBinI16(buffer, 16);
    int32_t nelements = readBinI32(buffer, 18);
    int32_t offset = readBinI32(buffer, 26);
    //int32_t handle = readBinI32(buffer, 30);
    //std::cout << filetype << '\t' << version << '\t' << name << '\t' << number << '\t' << etype << '\t' << esize << '\t' << nelements << '\t' << offset << '\t' << handle << std::endl;

    // Get all ABIF records
    std::vector<Abif> abi;
    for (int32_t i = 0; i < nelements; ++i) {
      int32_t ofs = i * esize + offset;
      std::vector<char> entry(buffer.begin()+ofs, buffer.begin()+ofs+esize);
      Abif ab;
      ab.name = readBinStr(entry, 0, 4);
      ab.number = readBinI32(entry, 4);
      ab.etype = readBinI16(entry, 8);
      ab.esize = readBinI16(entry, 10);
      ab.nelements = readBinI32(entry, 12);
      ab.dsize = readBinI32(entry, 16);
      ab.doffset = readBinI32(entry, 20);
      ab.key = ab.name + "." + _toString(ab.number);
      if (ab.name == "PCON") ab.etype = 1;
      abi.push_back(ab);
      //std::cout << ab.key << "\t" << ab.name << "\t" << ab.number << "\t" << ab.etype << "\t" << ab.esize << "\t" << ab.nelements << "\t" << ab.dsize << "\t" << ab.doffset << std::endl;
    }

    // Get what we need and dump the rest of this stupid format
    for(uint32_t i = 0; i < abi.size(); ++i) {
      int32_t ofs = i * esize + offset;
      int32_t ofsraw = ofs + 20;
      if (abi[i].dsize > 4) ofsraw = abi[i].doffset;
      std::vector<char> entry(buffer.begin()+ofsraw, buffer.begin()+ofsraw + abi[i].nelements*abi[i].esize + 1);
      if (abi[i].etype == 2) {
	if (abi[i].key == "PBAS.2") tr.basecalls1 = replaceNonDna(readBinStr(entry, 0, entry.size()));
	else if (abi[i].key == "P2BA.1") tr.basecalls2 = replaceNonDna(readBinStr(entry, 0, entry.size()));
	else if (abi[i].key == "FWO_.1") acgtOrder = readBinStr(entry, 0, entry.size());
      } else if (abi[i].etype == 4) {
	if (abi[i].key == "PLOC.2") {
	  for(int32_t k = 0; k < abi[i].nelements; ++k) {
	    tr.basecallpos.push_back(readBinI16(entry, k*2));
	  }
	} else if (abi[i].key == "DATA.9") {
	  for(int32_t k = 0; k < abi[i].nelements; ++k) {
	    trace[0].push_back(readBinI16(entry, k*2));
	  }
	} else if (abi[i].key == "DATA.10") {
	  for(int32_t k = 0; k < abi[i].nelements; ++k) {
	    trace[1].push_back(readBinI16(entry, k*2));
	  }
	} else if (abi[i].key == "DATA.11") {
	  for(int32_t k = 0; k < abi[i].nelements; ++k) {
	    trace[2].push_back(readBinI16(entry, k*2));
	  }
	} else if (abi[i].key == "DATA.12") {
	  for(int32_t k = 0; k < abi[i].nelements; ++k) {
	    trace[3].push_back(readBinI16(entry, k*2));
	  }
	}
      } else if (abi[i].etype == 1) {
	if (abi[i].key == "PCON.2") {
	  for(int32_t k = 0; k < abi[i].nelements; ++k) {
	    tr.qual.push_back(readBinUI8(entry, k));
	  }
	}
      }
    }
  }
  bfile.close();

  // Fix size of basecall vectors
  uint32_t minsize1 = tr.basecalls1.size();
  if (tr.basecalls2.size()) minsize1 = std::min(tr.basecalls1.size(), tr.basecalls2.size());
  uint32_t minsize2 = std::min(tr.qual.size(), tr.basecallpos.size());
  uint32_t minsize = std::min(minsize1, minsize2);
  tr.basecallpos.resize(minsize);
  tr.basecalls1.resize(minsize);
  tr.basecalls2.resize(minsize);
  tr.qual.resize(minsize);

  // Assign trace
  tr.traceACGT.resize(4, TMountains());
  for(uint32_t i = 0; i < acgtOrder.size(); ++i) {
    if (acgtOrder[i] == 'A') tr.traceACGT[0] = trace[i];
    else if (acgtOrder[i] == 'C') tr.traceACGT[1] = trace[i];
    else if (acgtOrder[i] == 'G') tr.traceACGT[2] = trace[i];
    else if (acgtOrder[i] == 'T') tr.traceACGT[3] = trace[i];
  }
  
  // Close input file
  return true;
}


inline void
basecall(Trace const& tr, BaseCalls& bc, float sigratio) {
  typedef Trace::TMountains TMountains;
  typedef Trace::TValue TValue;
  
  // Get peak regions
  std::vector<float> st;
  std::vector<float> ed;
  TValue oldVal = 0;
  TValue lastDiff = 0;
  for(uint32_t i = 0; i<tr.basecallpos.size(); ++i) {
    lastDiff = tr.basecallpos[i] - oldVal;
    st.push_back((float) tr.basecallpos[i] - 0.5 * (float) lastDiff);
    if (i > 0) ed.push_back((float) tr.basecallpos[i-1] + 0.5 * (float) lastDiff);
    oldVal = tr.basecallpos[i];
  }
  ed.push_back(tr.basecallpos[tr.basecallpos.size()-1] + 0.5 * lastDiff);

  // Call peaks
  std::vector<char> primary;
  std::vector<char> secondary;
  std::vector<char> consensus;
  for(uint32_t i = 0; i<st.size(); ++i) {
    TMountains pVal;
    TMountains pIdx;
    peak(tr.traceACGT, st[i], ed[i], pVal, pIdx);
    if ((pVal[0] == 0) && (pVal[1] == 0) && (pVal[2] == 0) && (pVal[3] == 0)) continue;
    TValue maxVal = 0;
    for(uint32_t k = 0; k<4; ++k)
      if (pVal[k] > maxVal) maxVal = pVal[k];
    std::vector<float> srat;
    for(uint32_t k = 0; k<4; ++k)
      srat.push_back((float) pVal[k] / (float) maxVal);
    float bestRat = sigratio;
    int32_t selACGT = -1;
    TValue selPos = 0;
    int32_t validBases = 0;
    for(uint32_t k = 0; k<4; ++k) {
      if (srat[k] >= sigratio) {
	++validBases;
	if (srat[k] >= bestRat) {
	  bestRat = srat[k];
	  selPos = pIdx[k];
	  selACGT = k;
	}
      }
    }
    //std::cout << bc.bcPos.size() << ',' << primary.size() << ',' << validBases << ',' << selPos << std::endl;
    bc.bcPos.push_back(selPos);
    if ((validBases == 4) || (selACGT == -1)) {
      primary.push_back('N');
      secondary.push_back('N');
      consensus.push_back('N');
    } else if (validBases > 1) {
      if (selACGT == 0) primary.push_back('A');
      else if (selACGT == 1) primary.push_back('C');
      else if (selACGT == 2) primary.push_back('G');
      else if (selACGT == 3) primary.push_back('T');
      TMountains leftover;
      for(int32_t k = 0; k<4; ++k) 
	if ((k != selACGT) && (srat[k] >= sigratio)) leftover.push_back(k);
      secondary.push_back(iupac(leftover));
      consensus.push_back('N');
    } else {
      if (selACGT == 0) {
	primary.push_back('A');
	secondary.push_back('A');
	consensus.push_back('A');
      } else if (selACGT == 1) {
	primary.push_back('C');
	secondary.push_back('C');
	consensus.push_back('C');
      } else if (selACGT == 2) {
	primary.push_back('G');
	secondary.push_back('G');
	consensus.push_back('G');
      } else if (selACGT == 3) {
	primary.push_back('T');
	secondary.push_back('T');
	consensus.push_back('T');
      }
    }
  }
  bc.primary = std::string(primary.begin(), primary.end());
  bc.secondary = std::string(secondary.begin(), secondary.end());
  bc.consensus = std::string(consensus.begin(), consensus.end());
}

inline void
traceTxtOut(std::string const& outfile, BaseCalls& bc, Trace const& tr) {
  typedef Trace::TValue TValue;
  uint32_t bcpos = 0;
  TValue idx = bc.bcPos[bcpos];
  std::ofstream rfile(outfile.c_str());
  rfile << "pos\tpeakA\tpeakC\tpeakG\tpeakT\tbasenum\tprimary\tsecondary\tconsensus\tqual" << std::endl;
  for(int32_t i = 0; i < (int32_t) tr.traceACGT[0].size(); ++i) {
    rfile << (i+1) << "\t";
    for(uint32_t k =0; k<4; ++k) rfile << tr.traceACGT[k][i] << "\t";
    if (idx == i) {
      rfile << (bcpos+1) << "\t";
      rfile << bc.primary[bcpos] << "\t" << bc.secondary[bcpos] << "\t" << bc.consensus[bcpos] << "\t" << (int32_t) tr.qual[bcpos] << "\t" << std::endl;
      if (bcpos < bc.bcPos.size() - 1) idx = bc.bcPos[++bcpos];
    } else rfile << "NA\tNA\tNA\tNA\tNA" << std::endl;
  }
}

template<typename TStream>
inline void
_traceJsonOut(TStream& rfile, BaseCalls& bc, Trace const& tr) {
  typedef Trace::TValue TValue;

  rfile << "\"pos\": [";
  for(uint32_t i = 0; i<tr.traceACGT[0].size(); ++i) {
    if (i!=0) rfile << ", ";
    rfile << (i+1);
  }
  rfile << "]," << std::endl;
  rfile << "\"peakA\": [";
  for(uint32_t i = 0; i<tr.traceACGT[0].size(); ++i) {
    if (i!=0) rfile << ", ";
    rfile << tr.traceACGT[0][i];
  }
  rfile << "]," << std::endl;
  rfile << "\"peakC\": [";
  for(uint32_t i = 0; i<tr.traceACGT[0].size(); ++i) {
    if (i!=0) rfile << ", ";
    rfile << tr.traceACGT[1][i];
  }
  rfile << "]," << std::endl;
  rfile << "\"peakG\": [";
  for(uint32_t i = 0; i<tr.traceACGT[0].size(); ++i) {
    if (i!=0) rfile << ", ";
    rfile << tr.traceACGT[2][i];
  }
  rfile << "]," << std::endl;
  rfile << "\"peakT\": [";
  for(uint32_t i = 0; i<tr.traceACGT[0].size(); ++i) {
    if (i!=0) rfile << ", ";
    rfile << tr.traceACGT[3][i];
  }
  rfile << "]," << std::endl;

  // Basecalls
  uint32_t bcpos = 0;
  TValue idx = bc.bcPos[0];
  rfile << "\"basecallPos\": [";
  for(int32_t i = 0; i < (int32_t) tr.traceACGT[0].size(); ++i) {
    if (idx == i) {
      if (i!=bc.bcPos[0]) rfile << ", ";
      rfile << (i+1);
      if (bcpos < bc.bcPos.size() - 1) idx = bc.bcPos[++bcpos];
    }
  }
  rfile << "]," << std::endl;
  bcpos = 0;
  idx = bc.bcPos[0];
  rfile << "\"basecalls\": {";
  for(int32_t i = 0; i < (int32_t) tr.traceACGT[0].size(); ++i) {
    if (idx == i) {
      if (i!=bc.bcPos[0]) rfile << ", ";
      rfile << "\"" << (i+1) << "\"" << ":" << "\"" << (bcpos+1) << ":" <<  bc.primary[bcpos];
      if (bc.primary[bcpos] != bc.secondary[bcpos]) rfile << "|" << bc.secondary[bcpos];
      rfile << "\"";
      if (bcpos < bc.bcPos.size() - 1) idx = bc.bcPos[++bcpos];
    }
  }
  rfile << "}" << std::endl;
}
   
inline void
traceJsonOut(std::string const& outfile, BaseCalls& bc, Trace const& tr) {
  // Output trace
  std::ofstream rfile(outfile.c_str());
  rfile << "{" << std::endl;
  _traceJsonOut(rfile, bc, tr);
  rfile << "}" << std::endl;
  rfile.close();  
}


}

#endif
