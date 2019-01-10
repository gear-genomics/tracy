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

#include "nlohmann/json.hpp"

using boost::asio::ip::tcp;

namespace tracy {

  inline int32_t
  variantsInRegion(std::string const& locus, std::string& respstr) {
    try {
      boost::asio::io_service io_service;

      // Connection
      std::string host("rest.ensembl.org");
      tcp::resolver resolver(io_service);
      tcp::resolver::query query(host, "http");
      tcp::resolver::iterator endpoint_iterator = resolver.resolve(query);
      tcp::socket socket(io_service);
      boost::asio::connect(socket, endpoint_iterator);
      
      // Request
      boost::asio::streambuf request;
      std::ostream request_stream(&request);

      // Debug
      request_stream << "GET " << "/overlap/region/human/" << locus << "?feature=variation; HTTP/1.0\r\n";
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


  inline int32_t
  parseJSON(std::string const& jsonString) {
    auto j = nlohmann::json::parse(jsonString);
    //std::cout << std::setw(4) << j << std::endl;

    // Collect annotated variants
    for(auto it = j.begin(); it != j.end(); ++it) {
      //std::cout << std::setw(4) << *it << std::endl;
      auto var = *it;
      if (var["alleles"].size() > 1) {
	std::string ref = var["alleles"][0];
	std::string alt = var["alleles"][1];
	std::string chr = var["seq_region_name"];
	std::string ct = var["consequence_type"];
	int32_t start = var["start"];
	int32_t end = var["end"];
	int32_t strand = var["strand"];
	std::string id = var["id"];
	std::cout << chr << ',' << start << ',' << ref << ',' << alt << std::endl;
      }
    }

    return 0;
  }
    

}

#endif
