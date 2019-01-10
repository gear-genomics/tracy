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

#ifndef VARIANTS_H
#define VARIANTS_H

#include <htslib/sam.h>
#include <htslib/vcf.h>

namespace tracy {

  struct Variant {
    int32_t pos;
    int32_t basenum;
    int32_t gt;
    std::string chr;
    std::string ref;
    std::string alt;
    std::string id;
    
    Variant(int32_t const p, int32_t const bc, int32_t const g, std::string const& c, std::string const& r, std::string const& a) : pos(p), basenum(bc), gt(g), chr(c), ref(r), alt(a), id(".") {}
  };

  template<typename TVar>
  struct SortVariant : public std::binary_function<TVar, TVar, bool>
  {
    inline bool operator()(TVar const& var1, TVar const& var2) {
      return ((var1.chr < var2.chr) || ((var1.chr == var2.chr) && (var1.pos < var2.pos)) || ((var1.chr == var2.chr) && (var1.pos == var2.pos) && (var1.basenum < var2.basenum)));
    }
  };
  
  
  inline void
  insertVariant(std::vector<Variant>& var, int32_t const pos, int32_t const bc, int32_t const gt, std::string const& chr, std::string const& ref, std::string const& alt) {
    // Search existing variant
    int32_t idx = -1;
    for(uint32_t i = 0; i < var.size(); ++i) {
      if ((var[i].pos == pos) && (var[i].chr == chr) && (var[i].ref == ref) && (var[i].alt == alt)) {
	idx = i;
	break;
      }
    }
    if (idx != -1) {
      // Update GT, homozygous variant
      var[idx].gt += 1;
    } else {
      // Insert variant
      var.push_back(Variant(pos, bc, gt, chr, ref, alt));
    }
  }

  
  template<typename TAlign>
  inline void
  callVariants(TAlign const& align, ReferenceSlice const& rs, std::vector<Variant>& var) {
    int32_t ri = rs.pos;
    
    // Get leading and trailing gaps
    int32_t viStart = 0;
    int32_t viEnd = 0;
    for(uint32_t j = 0; j < align.shape()[1]; ++j) {
      if (align[0][j] != '-') {
	if (viStart == 0) viStart = j;
	viEnd = j;
      }
      if ((align[1][j] != '-') && (viStart == 0)) ++ri;
    }

    // Call variants
    int32_t vi = 0;
    std::string del = "";
    int32_t delStart = 0;
    std::string ins = "";
    int32_t insStart = 0;
    char lastRefChar = '-';
    for(int32_t j = viStart; j <= viEnd; ++j) {
      if ((!del.empty()) && (align[0][j] != '-')) {
	// End of deletion
	insertVariant(var, delStart, vi, 1, rs.chr, del, std::string(1, del[0]));
	del = "";
      }
      if ((!ins.empty()) && (align[1][j] != '-')) {
	// End of insertion
	insertVariant(var, insStart, vi, 1, rs.chr, std::string(1, ins[0]), ins);
	ins = "";
      }
	
      if (align[0][j] != '-') ++vi;
      if (align[1][j] != '-') ++ri;
      if (align[0][j] != align[1][j]) {
	if ((align[0][j] != '-') && (align[1][j] != '-')) {
	  // SNV
	  insertVariant(var, ri, vi, 1, rs.chr, std::string(1, align[1][j]), std::string(1, align[0][j]));
	} else {
	  if (align[0][j] == '-') {
	    // Deletion
	    if (del.empty()) {
	      del.append(1, lastRefChar);
	      delStart = ri - 1;
	    }
	    del.append(1, align[1][j]);
	  } else {
	    // Insertion
	    if (ins.empty()) {
	      ins.append(1, lastRefChar);	      
	      insStart = ri;
	    }
	    ins.append(1, align[0][j]);
	  }	  
	}
      }
      if (align[1][j] != '-') lastRefChar = align[1][j];
    }

    
    //for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
    //std::cout << j << ':';
    //for(TAIndex i = 0; i < (TAIndex) align.shape()[0]; ++i) {
    //	std::cout << align[i][j];
    //}
    //std::cout << std::endl;
    //}
  }


  inline std::string
  variantType(std::string const& ref, std::string const& alt) {
    if ((ref.size() == 1) && (alt.size() == 1)) {
      return "SNV";
    } else {
      if (ref.size() > alt.size()) return "Deletion";
      else if (ref.size() < alt.size()) return "Insertion";
      else return "Complex";
    }
  }
  

  template<typename TConfig>
  inline void
  vcfOutput(TConfig const& c, BaseCalls const& bc, std::vector<Variant> const& var, ReferenceSlice const& rs) {
    // Output file name
    std::string outfile = c.outfile.string();
    if (c.format == "align") outfile = c.outfile.string() + ".bcf";
    
    // Output all structural variants
    htsFile *fp = hts_open(outfile.c_str(), "wb");
    bcf_hdr_t *hdr = bcf_hdr_init("w");

    // Print vcf header
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    boost::gregorian::date today = now.date();
    std::string datestr("##fileDate=");
    datestr += boost::gregorian::to_iso_string(today);
    bcf_hdr_append(hdr, datestr.c_str());
    bcf_hdr_append(hdr, "##FILTER=<ID=LowQual,Description=\"Low quality variant call.\">");
    bcf_hdr_append(hdr, "##INFO=<ID=ALIGNPOS,Number=1,Type=Integer,Description=\"Alignment position\">");
    bcf_hdr_append(hdr, "##INFO=<ID=BASEPOS,Number=1,Type=Integer,Description=\"Basecall position in trace\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SIGNALPOS,Number=1,Type=Integer,Description=\"Trace signal position\">");
    bcf_hdr_append(hdr, "##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Variant type\">");
    bcf_hdr_append(hdr, "##INFO=<ID=METHOD,Number=1,Type=String,Description=\"Type of approach used to detect variant\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");

    // Add reference
    std::string refloc("##reference=");
    refloc += c.genome.string();
    bcf_hdr_append(hdr, refloc.c_str());

    // Add contig info
    if (rs.filetype) {
      std::string refname("##contig=<ID=");
      refname += rs.chr + ",length=" + boost::lexical_cast<std::string>(rs.refslice.size()) + ">";
      bcf_hdr_append(hdr, refname.c_str());
    } else {
      faidx_t* fai = fai_load(c.genome.string().c_str());
      for(int32_t refIndex = 0; refIndex < faidx_nseq(fai); ++refIndex) {
	std::string seqname(faidx_iseq(fai, refIndex));
        uint32_t seqlen = faidx_seq_len(fai, seqname.c_str()) + 1;
	std::string refname("##contig=<ID=");
	refname += seqname + ",length=" + boost::lexical_cast<std::string>(seqlen) + ">";
	bcf_hdr_append(hdr, refname.c_str());
      }
      fai_destroy(fai);
    }
    
    // Add samples
    std::string sampleName = "sample";
    bcf_hdr_add_sample(hdr, sampleName.c_str());
    bcf_hdr_add_sample(hdr, NULL);
    bcf_hdr_write(fp, hdr);

    if (!var.empty()) {
      // Genotype arrays
      int32_t *gts = (int*) malloc(bcf_hdr_nsamples(hdr) * 2 * sizeof(int));
      int32_t *gqval = (int*) malloc(bcf_hdr_nsamples(hdr) * sizeof(int));
    
      // Iterate all structural variants
      bcf1_t *rec = bcf_init();
      for(uint32_t i = 0; i < var.size(); ++i) {

	// Output main vcf fields
	rec->rid = bcf_hdr_name2id(hdr, var[i].chr.c_str());
	rec->pos = var[i].pos;
	rec->qual = (int) bc.estQual[var[i].basenum];
	std::string id(var[i].id);
	bcf_update_id(hdr, rec, id.c_str());
	std::string alleles = var[i].ref + "," + var[i].alt;
	bcf_update_alleles_str(hdr, rec, alleles.c_str());
	int32_t tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");
	if (rec->qual < c.qualCut) tmpi = bcf_hdr_id2int(hdr, BCF_DT_ID, "LowQual");
	bcf_update_filter(hdr, rec, &tmpi, 1);
      
	// Add INFO fields
	std::string vartype = variantType(var[i].ref, var[i].alt);
	bcf_update_info_string(hdr, rec, "TYPE", vartype.c_str());
	std::string tracyVersion("EMBL.TRACYv");
	tracyVersion += tracyVersionNumber;
	bcf_update_info_string(hdr,rec, "METHOD", tracyVersion.c_str());
	tmpi = var[i].basenum;
	bcf_update_info_int32(hdr, rec, "ALIGNPOS", &tmpi, 1);
	tmpi = c.trimLeft + var[i].basenum;
	bcf_update_info_int32(hdr, rec, "BASEPOS", &tmpi, 1);
	tmpi = bc.bcPos[c.trimLeft + var[i].basenum];
	bcf_update_info_int32(hdr, rec, "SIGNALPOS", &tmpi, 1);
	
	// Add genotyping information
	if (var[i].gt == 0) {
	  gts[0] = bcf_gt_unphased(0);
	  gts[1] = bcf_gt_unphased(0);
	} else if (var[i].gt == 1) {
	  gts[0] = bcf_gt_unphased(0);
	  gts[1] = bcf_gt_unphased(1);
	} else if (var[i].gt == 2) {
	  gts[0] = bcf_gt_unphased(1);
	  gts[1] = bcf_gt_unphased(1);
	} else {
	  gts[0] = bcf_gt_missing;
	  gts[1] = bcf_gt_missing;
	}
	gqval[0] = (int) bc.estQual[var[i].basenum];
	bcf_update_genotypes(hdr, rec, gts, bcf_hdr_nsamples(hdr) * 2);
	bcf_update_format_int32(hdr, rec, "GQ", gqval, bcf_hdr_nsamples(hdr));

	// Write record
	bcf_write1(fp, hdr, rec);
	bcf_clear1(rec);
      }
      bcf_destroy1(rec);
    
      // Clean-up
      free(gts);
      free(gqval);
    }
    // Close VCF file
    bcf_hdr_destroy(hdr);
    hts_close(fp);
  
    // Build index
    bcf_index_build(outfile.c_str(), 14);
  }


}

#endif
