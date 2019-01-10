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

#ifndef WEB_H
#define WEB_H

#include <boost/asio.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <json.hpp>

using boost::asio::ip::tcp;

namespace tracy {


  struct KnownVariation {
    typedef std::vector<std::string> TAlleles;
    int32_t pos;
    std::string id;
    std::string chr;
    TAlleles alleles;

    KnownVariation(int32_t const p, std::string const& name, std::string const& c, std::vector<std::string> const& a) : pos(p), id(name), chr(c), alleles(a) {}
  };

  std::string speciesArray[] = {"ciona_savignyi", "erinaceus_europaeus", "mastacembelus_armatus", "cavia_porcellus", "homo_sapiens", "pelodiscus_sinensis", "meleagris_gallopavo", "chlorocebus_sabaeus", "sphenodon_punctatus", "saimiri_boliviensis_boliviensis", "carlito_syrichta", "poecilia_mexicana", "chrysemys_picta_bellii", "poecilia_formosa", "mus_musculus_lpj", "cercocebus_atys", "tupaia_belangeri", "saccharomyces_cerevisiae", "oryzias_latipes", "tursiops_truncatus", "papio_anubis", "propithecus_coquereli", "myotis_lucifugus", "pan_paniscus", "mus_musculus_129s1svimj", "gorilla_gorilla", "amphilophus_citrinellus", "maylandia_zebra", "rhinopithecus_roxellana", "bos_taurus", "oryctolagus_cuniculus", "poecilia_reticulata", "otolemur_garnettii", "mus_musculus", "mus_musculus_casteij", "mesocricetus_auratus", "cebus_capucinus", "drosophila_melanogaster", "mola_mola", "sarcophilus_harrisii", "gasterosteus_aculeatus", "cricetulus_griseus_chok1gshd", "loxodonta_africana", "chinchilla_lanigera", "felis_catus", "mus_musculus_cbaj", "sorex_araneus", "taeniopygia_guttata", "heterocephalus_glaber_female", "mus_musculus_balbcj", "nomascus_leucogenys", "rhinopithecus_bieti", "mus_caroli", "mus_musculus_nodshiltj", "haplochromis_burtoni", "pongo_abelii", "heterocephalus_glaber_male", "danio_rerio", "caenorhabditis_elegans", "seriola_lalandi_dorsalis", "eptatretus_burgeri", "tetraodon_nigroviridis", "oryzias_melastigma", "fundulus_heteroclitus", "canis_lupus_dingo", "ficedula_albicollis", "xiphophorus_couchianus", "mustela_putorius_furo", "mus_musculus_dba2j", "acanthochromis_polyacanthus", "hippocampus_comes", "mus_spretus", "pygocentrus_nattereri", "amphiprion_ocellaris", "mus_musculus_akrj", "takifugu_rubripes", "procavia_capensis", "oreochromis_niloticus", "latimeria_chalumnae", "astyanax_mexicanus", "labrus_bergylta", "aotus_nancymaae", "seriola_dumerili", "stegastes_partitus", "ovis_aries", "cricetulus_griseus_crigri", "canis_familiaris", "mus_musculus_fvbnj", "gambusia_affinis", "amphiprion_percula", "periophthalmus_magnuspinnatus", "mus_musculus_aj", "vulpes_vulpes", "equus_asinus_asinus", "octodon_degus", "callithrix_jacchus", "mandrillus_leucophaeus", "ciona_intestinalis", "mus_pahari", "anabas_testudineus", "rattus_norvegicus", "notamacropus_eugenii", "monodelphis_domestica", "equus_caballus", "gopherus_agassizii", "panthera_pardus", "petromyzon_marinus", "scleropages_formosus", "microcebus_murinus", "pundamilia_nyererei", "anas_platyrhynchos", "astatotilapia_calliptera", "ictalurus_punctatus", "capra_hircus", "dipodomys_ordii", "ursus_maritimus", "macaca_mulatta", "pan_troglodytes", "poecilia_latipinna", "ursus_americanus", "macaca_fascicularis", "cyprinodon_variegatus", "peromyscus_maniculatus_bairdii", "choloepus_hoffmanni", "xiphophorus_maculatus", "nannospalax_galili", "panthera_tigris_altaica", "anolis_carolinensis", "mus_musculus_wsbeij", "microtus_ochrogaster", "jaculus_jaculus", "ornithorhynchus_anatinus", "phascolarctos_cinereus", "scophthalmus_maximus", "sus_scrofa", "mus_musculus_pwkphj", "ailuropoda_melanoleuca", "fukomys_damarensis", "esox_lucius", "xenopus_tropicalis", "cynoglossus_semilaevis", "dasypus_novemcinctus", "gallus_gallus", "oryzias_latipes_hni", "paramormyrops_kingsleyae", "vicugna_pacos", "mus_musculus_nzohlltj", "colobus_angolensis_palliatus", "monopterus_albus", "kryptolebias_marmoratus", "oryzias_latipes_hsok", "mus_musculus_c3hhej", "pteropus_vampyrus", "neolamprologus_brichardi", "cavia_aperea", "gadus_morhua", "lepisosteus_oculatus", "macaca_nemestrina", "echinops_telfairi", "ochotona_princeps", "ictidomys_tridecemlineatus", "mus_musculus_c57bl6nj"};
  
  inline bool
  speciesExist(std::string const& tiger) {
    //curl -o species -X GET --header 'Accept: application/json' http://rest.ensembl.org/info/species
    typedef std::set<std::string> TSpecies;
    TSpecies species(speciesArray, speciesArray + sizeof(speciesArray) / sizeof(speciesArray[0]));
    if (species.find(tiger) != species.end()) return true;
    else return false;
  }


  template<typename TConfig>
  inline int32_t
    variantsInRegion(TConfig const& c, std::string const& locus, std::string& respstr) {
    try {
      boost::asio::io_service io_service;

      // Connection
      std::string host(c.host);
      tcp::resolver resolver(io_service);
      tcp::resolver::query query(host, "http");
      tcp::resolver::iterator endpoint_iterator = resolver.resolve(query);
      tcp::socket socket(io_service);
      boost::asio::connect(socket, endpoint_iterator);
      
      // Request
      boost::asio::streambuf request;
      std::ostream request_stream(&request);

      // Debug
      request_stream << "GET " << "/overlap/region/" << c.annotate << "/" << locus << "?feature=variation; HTTP/1.0\r\n";
      request_stream << "Host: " << host << "\r\n";
      request_stream << "Accept: application/json\r\n";
      request_stream << "Connection: close\r\n\r\n";
      boost::asio::write(socket, request);
      
      // Respose
      boost::asio::streambuf response;
      boost::asio::read_until(socket, response, "\r\n");
      
      // Check response
      std::istream response_stream(&response);
      std::string http_version;
      response_stream >> http_version;
      unsigned int status_code;
      response_stream >> status_code;
      std::string status_message;
      std::getline(response_stream, status_message);
      if ((!response_stream) || (http_version.substr(0, 5) != "HTTP/")) {
	std::cerr << "Invalid response" << std::endl;
	return 1;
      }
      if (status_code != 200) {
	std::cerr << "Response returned with status code " << status_code << "\n";
	return 1;
      }
      
      // Response header
      boost::system::error_code error;
      boost::asio::read_until(socket, response, "\r\n\r\n");
      std::string header;
      while(std::getline(response_stream, header) && header != "\r") {
	// Debug
	//std::cerr << header << std::endl;
      }
      
      // Content
      if (response.size() > 0) {
	std::ostringstream rstream;
	rstream << &response;
	respstr += rstream.str();
      }
      
      // Until EOF
      while(boost::asio::read(socket, response, boost::asio::transfer_at_least(1), error)) {
	if (response.size() > 0) {
	  std::ostringstream rstream;
	  rstream << &response;
	  respstr += rstream.str();
	}
      }
      if (error != boost::asio::error::eof) {
	throw boost::system::system_error(error);
      }
    } catch (std::exception& e) {
      std::cerr << "Exception: " << e.what() << std::endl;
    }
      
    return 0;
  }

  inline void
  annotateVariants(std::vector<KnownVariation> const& kv, std::vector<Variant>& var) {
    for(uint32_t i = 0; i < var.size(); ++i) {
      for(uint32_t j = 0; j < kv.size(); ++j) {
	if (var[i].pos == kv[j].pos) {
	  if (var[i].chr == kv[j].chr) {
	    // Debug
	    //std::cout << var[i].ref << ',' << var[i].alt << std::endl;
	    //for(uint32_t k = 0; k < kv[j].alleles.size(); ++k) {
	    //std::cout << kv[j].alleles[k] << ',';
	    //}
	    //std::cout << std::endl;
	    // Match alleles
	    if (var[i].ref == kv[j].alleles[0]) {
	      for(uint32_t k = 1; k < kv[j].alleles.size(); ++k) {
		if (var[i].alt == kv[j].alleles[k]) {
		  var[i].id = kv[j].id;
		  break;
		}
	      }
	    }	  
	  }
	}
      }
    }
  }

  inline int32_t
  parseKnownVariants(std::string const& jsonString, std::vector<KnownVariation>& kv) {
    auto j = nlohmann::json::parse(jsonString);
    //std::cout << std::setw(4) << j << std::endl;

    // Collect annotated variants
    for(auto it = j.begin(); it != j.end(); ++it) {
      //std::cout << std::setw(4) << *it << std::endl;
      auto var = *it;
      if (var["alleles"].size() > 1) {
	int32_t strand = var["strand"];
	// Fwd strand
	if (strand == 1) {
	  int32_t start = var["start"];
	  int32_t end = var["end"];
	  // Only SNPs
	  if (start == end) {
	    std::string chr = var["seq_region_name"];
	    std::string id = var["id"];
	    std::vector<std::string> alleles;
	    for(uint32_t i = 0; i <  var["alleles"].size(); ++i) {
	      if (var["alleles"][i] != "-") {
		std::string al = var["alleles"][i];
		if (al.size() == 1) alleles.push_back(al);
	      }
	    }
	    if (alleles.size() > 1) kv.push_back(KnownVariation(start, id, chr, alleles));
	  }
	}
      }
    }
    return kv.size();
  }
    

}

#endif
