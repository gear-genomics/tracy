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
      boost::asio::read_until(socket, response, "\r\n\r\n");
      std::string header;
      while(std::getline(response_stream, header) && header != "\r") {
	// Debug
	//std::cerr << header << std::endl;
      }
      
      // Content
      if (response.size() > 0) {
	std::istream rstream(&response);
	std::string resplocal;
	rstream >> resplocal;
	respstr += resplocal;
      }

      // Until EOF
      boost::system::error_code error;
      while(boost::asio::read(socket, response, boost::asio::transfer_at_least(1), error)) {
	std::istream rstream(&response);
	std::string resplocal;
	rstream >> resplocal;
	respstr += resplocal;
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
  parseJSON(std::string const& json) {
    //Debug
    //std::cerr << json << std::endl;
    
    std::stringstream ss;
    ss << json;
    boost::property_tree::ptree root;
    boost::property_tree::read_json(ss, root);
    return 0;
  }
    

}

#endif
