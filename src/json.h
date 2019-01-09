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

#ifndef JSON_H
#define JSON_H

#include <boost/progress.hpp>
#include "abif.h"
#include "fmindex.h"
#include "version.h"
#include "variants.h"

namespace tracy
{
  
  #ifndef EMPTY_TRACE_SIGNAL
  #define EMPTY_TRACE_SIGNAL -99
  #endif  


  template<typename TStream, typename TConfig>
  inline void
  _metaOut(TStream& rfile, TConfig const& c) {
    rfile << "\"meta\": {";
    rfile << "\"program\": \"tracy\", ";
    rfile << "\"version\": \"" << tracyVersionNumber << "\", ";
    rfile << "\"arguments\": {";
    rfile << "\"trimLeft\": " << c.trimLeft << ", ";
    rfile << "\"trimRight\": " << c.trimRight << ", ";
    rfile << "\"pratio\": " << c.pratio << ", ";
    rfile << "\"genome\": \"" << c.genome.string() << "\", ";
    rfile << "\"input\": \"" << c.ab.string() << "\"";
    rfile << "}}," << std::endl;
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
	if (bc.primary[bcpos] != bc.secondary[bcpos]) rfile << "|" << expandIUPAC(bc.secondary[bcpos]);
	rfile << "\"";
	if (bcpos < bc.bcPos.size() - 1) idx = bc.bcPos[++bcpos];
      }
    }
    rfile << "}," << std::endl;
    rfile << "\"primarySeq\": \"" << bc.primary << "\"," << std::endl;
    rfile << "\"secondarySeq\": \"" << bc.secondary << "\"" << std::endl;
  }
  
  inline void
  traceJsonOut(std::string const& outfile, BaseCalls& bc, Trace const& tr) {
    // Output trace
    std::ofstream rfile(outfile.c_str());
    rfile << "{" << std::endl;
    _traceJsonOut(rfile, bc, tr);
    rfile << std::endl;
    rfile << "}" << std::endl;
    rfile.close();  
  }

  
  template<typename TAlign>
  inline void
  traceAlignJsonOut(std::string const& outfile, BaseCalls& bc, Trace const& tr, ReferenceSlice const& rs, TAlign const& align) {
    typedef Trace::TValue TValue;

    
    // Output trace
    std::ofstream rfile(outfile.c_str());
    rfile << "{" << std::endl;
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
    uint32_t gaplessbcpos = 0;
    rfile << "\"basecalls\": {";
    for(int32_t i = 0; i < (int32_t) tr.traceACGT[0].size(); ++i) {
      if (idx == i) {
	if (i!=bc.bcPos[0]) rfile << ", ";
	if (bc.primary[bcpos] != '-') {
	  rfile << "\"" << (i+1) << "\"" << ":" << "\"" << (++gaplessbcpos) << ":" <<  bc.primary[bcpos];
	  if (bc.primary[bcpos] != bc.secondary[bcpos]) rfile << "|" << bc.secondary[bcpos];
	  rfile << "\"";
	} else rfile << "\"" << (i+1) << "\"" << ":" << "\"-\"";
	if (bcpos < bc.bcPos.size() - 1) idx = bc.bcPos[++bcpos];
      }
    }
    rfile << "}," << std::endl;
    rfile << "\"refchr\": \"" << rs.chr << "\"," << std::endl;
    rfile << "\"refpos\": " << (rs.pos + 1) << "," << std::endl;
    rfile << "\"altalign\": \"";
    for(uint32_t j = 0; j<align.shape()[1]; ++j) rfile << align[0][j];
    rfile << "\"," << std::endl;
    rfile << "\"refalign\": \"";
    for(uint32_t j = 0; j<align.shape()[1]; ++j) rfile << align[1][j];
    rfile << "\"," << std::endl;
    rfile << "\"forward\": " << rs.forward << std::endl;
    rfile << "}" << std::endl;
    rfile.close();  
  }


  template<typename TConfig, typename TAlign, typename TDecomposition>
  inline void
  traceAlleleAlignJsonOut(TConfig const& c, BaseCalls& bc, Trace const& tr, std::vector<Variant> const& var, ReferenceSlice const& rs1, ReferenceSlice const& rs2, ReferenceSlice const&, TAlign const& align1, TAlign const& align2, TAlign const& align3, TDecomposition const& dcp, int32_t const a1Score, int32_t const a2Score, int32_t const a3Score, TraceBreakpoint const& bp, std::pair<double, double> const& a1a2) {
    // Output file name
    std::string outfile = c.outfile.string();
    if (c.format == "align") outfile = c.outfile.string() + ".json";    
    
    // Output trace
    std::ofstream rfile(outfile.c_str());
    rfile << "{" << std::endl;

    // Meta output
    _metaOut(rfile, c);
    
    // Trace Output
    _traceJsonOut(rfile, bc, tr);
    rfile << "," << std::endl;

    // Provide x-window
    int32_t lb = bc.bcPos[c.trimLeft + bp.breakpoint];
    if (lb < 300) lb = 1;
    else lb -= 300;
    int32_t ub = bc.bcPos[c.trimLeft + bp.breakpoint];
    if (ub + 300 < bc.bcPos[bc.bcPos.size() - 1]) ub += 300;
    else ub = bc.bcPos[bc.bcPos.size() - 1];
    rfile << "\"chartConfig\": { \"x\": { \"axis\": { \"range\": [" << lb << ", " << ub << "] }}}," << std::endl;
    
    // Allele1
    rfile << "\"ref1chr\": \"" << rs1.chr << "\"," << std::endl;
    rfile << "\"ref1pos\": " << (rs1.pos + 1) << "," << std::endl;
    rfile << "\"alt1align\": \"";
    for(uint32_t j = 0; j<align1.shape()[1]; ++j) rfile << align1[0][j];
    rfile << "\"," << std::endl;
    rfile << "\"ref1align\": \"";
    for(uint32_t j = 0; j<align1.shape()[1]; ++j) rfile << align1[1][j];
    rfile << "\"," << std::endl;
    rfile << "\"ref1forward\": " << rs1.forward << "," << std::endl;
    rfile << "\"align1score\": " << a1Score << "," << std::endl;

    // Allele 2
    rfile << "\"ref2chr\": \"" << rs2.chr << "\"," << std::endl;
    rfile << "\"ref2pos\": " << (rs2.pos + 1) << "," << std::endl;
    rfile << "\"alt2align\": \"";
    for(uint32_t j = 0; j<align2.shape()[1]; ++j) rfile << align2[0][j];
    rfile << "\"," << std::endl;
    rfile << "\"ref2align\": \"";
    for(uint32_t j = 0; j<align2.shape()[1]; ++j) rfile << align2[1][j];
    rfile << "\"," << std::endl;
    rfile << "\"ref2forward\": " << rs2.forward << "," << std::endl;
    rfile << "\"align2score\": " << a2Score << "," << std::endl;

    // Alt1 vs. Alt2
    rfile << "\"allele1fraction\": " << a1a2.first << "," << std::endl;
    rfile << "\"allele1align\": \"";
    for(uint32_t j = 0; j<align3.shape()[1]; ++j) rfile << align3[0][j];
    rfile << "\"," << std::endl;
    rfile << "\"allele2fraction\": " << a1a2.second << "," << std::endl;
    rfile << "\"allele2align\": \"";
    for(uint32_t j = 0; j<align3.shape()[1]; ++j) rfile << align3[1][j];
    rfile << "\"," << std::endl;
    rfile << "\"align3score\": " << a3Score << "," << std::endl;

    // Breakpoint
    rfile << "\"hetindel\": " << bp.indelshift << "," << std::endl;
    rfile << "\"hetindelbp\": " << c.trimLeft + bp.breakpoint << "," << std::endl;
    
    // Decomposition
    rfile << "\"decomposition\": " << "{" << std::endl;
    rfile << "\"x\": [";
    for(uint32_t i = 0; i < dcp.size(); ++i) {
      if (i!=0) rfile << ", ";
      rfile << dcp[i].first;
    }
    rfile << "]," << std::endl;
    rfile << "\"y\": [";
    for(uint32_t i = 0; i < dcp.size(); ++i) {
      if (i!=0) rfile << ", ";
      rfile << dcp[i].second;
    }
    rfile << "]" << std::endl;
    if (var.empty()) rfile << "}" << std::endl;
    else rfile << "}," << std::endl;

    // Variants
    if (!var.empty()) {
      rfile << "\"variants\": {" <<  std::endl;
      rfile << "\"columns\": [";
      rfile << "\"chr\", \"pos\", \"id\", \"ref\", \"alt\", \"qual\", \"filter\", \"type\", \"basenum\", \"genotype\"" << std::endl;
      rfile << "]," << std::endl;
      rfile << "\"rows\": [" << std::endl;
      for(uint32_t i = 0; i < var.size(); ++i) {
	if (i > 0) rfile << "," << std::endl;
	rfile << "[";
	rfile << "\"" << var[i].chr << "\", ";
	rfile << var[i].pos << ", ";
	rfile << "\".\", ";
	rfile << "\"" << var[i].ref << "\", ";
	rfile << "\"" << var[i].alt << "\", ";
	rfile << (int32_t) bc.estQual[var[i].basenum] << ", ";
	if ((int32_t) bc.estQual[var[i].basenum] < c.qualCut) rfile << "\"LowQual\", ";
	else rfile << "\"PASS\", ";
	rfile << "\"" << variantType(var[i].ref, var[i].alt) << "\", ";
	rfile << var[i].basenum << ", ";
	if (var[i].gt == 0) rfile << "\"hom. REF\"";
	else if (var[i].gt == 1) rfile << "\"het.\"";
	else if (var[i].gt == 2) rfile << "\"hom. ALT\"";
	else rfile << "\"missing\"";
	rfile << "]";
      }
      rfile << "]" << std::endl;
      rfile << "}" << std::endl;
    }
    
    // Close
    rfile << "}" << std::endl;
    rfile.close();
  }
  
  template<typename TAlign>
  inline void
  alignmentTracePadding(TAlign const& align, Trace const& tr, BaseCalls const& bc, Trace& ntr, BaseCalls& nbc) {
    typedef Trace::TMountains TMountains; 
    typedef Trace::TValue TValue;
    typedef std::vector<uint32_t> TVecVal;
    TVecVal insPos;
    TVecVal insSize; 

    // Calculate average trace basecall offset
    uint32_t step = 6;
    if (bc.bcPos.size()>1) {
      double avg = 0;
      for(uint32_t i = 1; i < bc.bcPos.size(); ++i) avg += (bc.bcPos[i] - bc.bcPos[i-1]);
      avg /= (bc.bcPos.size() - 1);
      step = (uint32_t) avg;
    }

    // Find the gap positions 
    uint32_t pos = 0;
    bool ingap = false;
    uint32_t gapsize = 0;    
    for(uint32_t j = 0; j<align.shape()[1]; ++j) {
      if (align[0][j] == '-') {
	if (ingap) ++gapsize;
	else {
	  gapsize = 1;
	  ingap = true;
	}
      } else {
	if (ingap) {
	  ingap = false;
	  if (pos) {
	    uint32_t insertPos = (uint32_t) (( bc.bcPos[pos - 1] + bc.bcPos[pos] ) / 2.0);
	    insPos.push_back(insertPos);
	    insSize.push_back(gapsize);
	  } else {
	    insPos.push_back(0);
	    insSize.push_back(gapsize);
	  }
	}
	++pos;
      }
    }
    // Trailing gaps
    if (ingap) {
      uint32_t insertPos = tr.traceACGT[0].size();
      insPos.push_back(insertPos);
      insSize.push_back(gapsize);
    }
    //for(uint32_t i = 0; i<insPos.size(); ++i) std::cerr << i << ',' << insPos[i] << ',' << insSize[i] << std::endl;
    

    // Rewrite the arrays
    uint32_t bcpos = 0;
    TValue idx = bc.bcPos[0];
    ntr.traceACGT.resize(4, TMountains());
    uint32_t offset = 0;
    TValue tracePos = 0;
    uint32_t inspos = 0;
    TValue insIdx = -1;
    if (!insPos.empty()) insIdx = insPos[0];
    for(; tracePos < (TValue) tr.traceACGT[0].size(); ++tracePos) {
      ntr.traceACGT[0].push_back(tr.traceACGT[0][tracePos]);        
      ntr.traceACGT[1].push_back(tr.traceACGT[1][tracePos]);
      ntr.traceACGT[2].push_back(tr.traceACGT[2][tracePos]);
      ntr.traceACGT[3].push_back(tr.traceACGT[3][tracePos]);
      if (insIdx == tracePos) {
	for(uint32_t k = 0; k < insSize[inspos]; ++k) {
	  nbc.bcPos.push_back(tracePos + offset + (uint32_t) (step / 2.0));
	  nbc.primary.push_back('-');
	  nbc.secondary.push_back('-');
	  nbc.consensus.push_back('-');
	  for(uint32_t n = 0; n < step; ++n, ++offset) {
	    ntr.traceACGT[0].push_back(EMPTY_TRACE_SIGNAL);
	    ntr.traceACGT[1].push_back(EMPTY_TRACE_SIGNAL);
	    ntr.traceACGT[2].push_back(EMPTY_TRACE_SIGNAL);
	    ntr.traceACGT[3].push_back(EMPTY_TRACE_SIGNAL);
	  }
        }
	if (inspos < insPos.size() - 1) insIdx = insPos[++inspos];
      }
      if (idx == tracePos) {
	nbc.bcPos.push_back(idx + offset);
	nbc.primary.push_back(bc.primary[bcpos]);
	nbc.secondary.push_back(bc.secondary[bcpos]);
	nbc.consensus.push_back(bc.consensus[bcpos]);
	if (bcpos < bc.bcPos.size() - 1) idx = bc.bcPos[++bcpos];
      }
    }
    // Trailing gaps
    if (insIdx == tracePos) {
      for(uint32_t k = 0; k < insSize[inspos]; ++k) {
	nbc.bcPos.push_back(tracePos + offset + (uint32_t) (step / 2.0));
	nbc.primary.push_back('-');
	nbc.secondary.push_back('-');
	nbc.consensus.push_back('-');
	for(uint32_t n = 0; n < step; ++n, ++offset) {
	  ntr.traceACGT[0].push_back(EMPTY_TRACE_SIGNAL);
	  ntr.traceACGT[1].push_back(EMPTY_TRACE_SIGNAL);
	  ntr.traceACGT[2].push_back(EMPTY_TRACE_SIGNAL);
	  ntr.traceACGT[3].push_back(EMPTY_TRACE_SIGNAL);
	}
      }
    }
  }

}

#endif
