/*-----------------------------------------------------------------------------*\
                  ___   _____   _____   _____     ___   
                 / ___\/\  __`\/\  __`\/\  __`\  / __`\ 
                /\ \__/\ \ \_\ \ \ \_\ \ \ \_\ \/\ \_\ \
                \ \____\\ \  __/\ \  __/\ \  __/\ \____/
                 \/____/ \ \ \/  \ \ \/  \ \ \/  \/___/ 
                          \ \_\   \ \_\   \ \_\         
                           \/_/    \/_/    \/_/         

         A Compilation for Fluid-Particle Data Post PrOcessing

Copyright (C): 2014 DCS Computing GmbH (www.dcs-computing.com), Linz, Austria
               2014 Graz University of Technology (ippt.tugraz.at), Graz, Austria
---------------------------------------------------------------------------------
License
    CPPPO is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with CPPPO. If not, see <http://www.gnu.org/licenses/lgpl.html>.

	This code is designed for on-the-fly post processing of fluid-particle
	data (e.g., of velocity, pressure, concentration, temperature field).

	Parts of the code were developed in the frame of the NanoSim project funded
	by the European Commission through FP7 Grant agreement no. 604656.
\*-----------------------------------------------------------------------------*/

#ifndef C3PO_OUTPUT_H
#define C3PO_OUTPUT_H

#include "stdio.h"
#include "c3po_base.h"
#include "c3po_base_interface.h"
#include "qjson_includes.h"


namespace C3PO_NS
{




class Output : public c3poBase, public c3poBaseInterface
{
    public:

      Output(c3po *ptr);
      ~Output();

      void write_screen_one(const char *message) const;
      void write_screen_all(const char *message) const;
      void write_screen_one_int(const int number, const char *message) const;
      void write_screen_one_int2(const int number, const char *message,const char *message2) const;
      void write_screen_all_int(const int number, const char *message) const;
       
      void write_log_one(const char *message) const;
      void write_log_all(const char *message) const;
      void write_time_one(double *message);
     
      static void createQJsonArrays(std::string,std::string,
                                    std::vector<std::string>, std::vector<double*>, 
                                    int datanum = -1, bool overwrite=true, 
                                    std::vector<int> * datanumvec_ = NULL
                                   );
                                   
      static void createQJsonArrays(std::string,std::string,
                                    std::vector<std::string>, std::vector<int*>, 
                                    int datanum = -1, bool overwrite=true, 
                                    std::vector<int> * datanumvec_ = NULL
                                   );


      static void createQJsonArrays(std::string,std::string,
                                    std::vector<std::string>,std::vector<double**>,
                                    int,int,int,bool
                                   );
 
      void generateDir(std::string dir, bool &check) const;

      //void write_screen_log_one(const char *message);
      //void write_screen_log_all(const char *message);

    private:

      FILE *screen_;                  // screen output
      FILE *logfile_;                 // logfile
      char *version_;                      // basic identifier of code version

};

} //end c3po_NS

#endif
