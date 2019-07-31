#ifndef VERSION_H
#define VERSION_H

namespace tracy
{


  std::string tracyVersionNumber = "0.5.3";

  inline 
    void printTitle(std::string const& title) 
    {
      std::cout << "**********************************************************************" << std::endl;
      std::cout << "Program: Tracy" << std::endl;
      std::cout << "This is free software, and you are welcome to redistribute it under" << std::endl;
      std::cout << "certain conditions (BSD License); for license details use '-l'." << std::endl;
      std::cout << "This program comes with ABSOLUTELY NO WARRANTY; for details use '-w'." << std::endl;
      std::cout <<  std::endl;
      std::cout <<  title << " (Version: " << tracyVersionNumber << ")" << std::endl;
      std::cout << "Contact: Tobias Rausch (rausch@embl.de)" << std::endl;
      std::cout << "**********************************************************************" << std::endl;
      std::cout << std::endl;
    }

  inline
    void displayWarranty()
    {
      std::cout << "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND" << std::endl;
      std::cout << "ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED" << std::endl;
      std::cout << "WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE" << std::endl;
      std::cout << "DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY" << std::endl;
      std::cout << "DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES" << std::endl;
      std::cout << "(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;" << std::endl;
      std::cout << "LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND" << std::endl;
      std::cout << "ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT" << std::endl;
      std::cout << "(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS" << std::endl;
      std::cout << "SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE." << std::endl;
      std::cout << std::endl;
    }
    
  inline void
  bsd() {
    std::cout << "Copyright (c) 2019, the Tracy Project Authors." << std::endl;
    std::cout << "All rights reserved. Please see the AUTHORS file for details." << std::endl;
    std::cout << std::endl;
    std::cout << "Redistribution and use in source and binary forms, with or without" << std::endl;
    std::cout << "modification, are permitted provided that the following conditions are met:" << std::endl;
    std::cout << "    * Redistributions of source code must retain the above copyright" << std::endl;
    std::cout << "      notice, this list of conditions and the following disclaimer." << std::endl;
    std::cout << "    * Redistributions in binary form must reproduce the above copyright" << std::endl;
    std::cout << "      notice, this list of conditions and the following disclaimer in the" << std::endl;
    std::cout << "      documentation and/or other materials provided with the distribution." << std::endl;
    std::cout << "    * Neither the name of the copyright holder nor the" << std::endl;
    std::cout << "      names of its contributors may be used to endorse or promote products" << std::endl;
    std::cout << "      derived from this software without specific prior written permission." << std::endl;
    std::cout << std::endl;
    std::cout << "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND" << std::endl;
    std::cout << "ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED" << std::endl;
    std::cout << "WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE" << std::endl;
    std::cout << "DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY" << std::endl;
    std::cout << "DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES" << std::endl;
    std::cout << "(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;" << std::endl;
    std::cout << "LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND" << std::endl;
    std::cout << "ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT" << std::endl;
    std::cout << "(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS" << std::endl;
    std::cout << "SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE." << std::endl;
    std::cout << std::endl;
  }

}

#endif
