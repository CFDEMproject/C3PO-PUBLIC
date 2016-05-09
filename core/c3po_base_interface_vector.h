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

#ifndef C3PO_BASE_INTERFACE_VECTOR_H
#define C3PO_BASE_INTERFACE_VECTOR_H

#include "c3po_base_interface.h"
#include "c3po.h"
#include "mpi.h"
#include <string>
#include <vector>
#include <map>


namespace C3PO_NS
{

class c3poBaseInterfaceVector
{
  friend class c3po;

  public:

    void init() const
    {
        // HOW C3PO WORKS

        // -1 make settings (any commands in input script)
        // -2 init (this function, called before a run is executed)
        // -3 run (via run command), including phyiscs, coupling, output
        // - can repeat 1-3 as often as needed

        // read all data; either from JSON files or restart files
        // performed on proc 0
        for (unsigned i=0; i < list_.size(); i++)
            list_[i]->read();

        // restart data from previous run if applicable
        // performed on proc 0
        for (unsigned i=0; i < list_.size(); i++)
            list_[i]->restart();

        // scatter all properties/settings so all MPI procs have it
        for (unsigned i=0; i < list_.size(); i++)
            list_[i]->scatter();

        // parallellize - each proc just keeps data he needs
        for (unsigned i=0; i < list_.size(); i++)
            list_[i]->parallelize();

        // do model/class specific init after everything is scattered
        for (unsigned i=0; i < list_.size(); i++)
            list_[i]->init();
    }

    void parse_command(std::string map_string,int narg,char const* const* arg) const
    {
        if (map_.find(map_string) != map_.end())
        {
            map_[map_string]->parse_command(narg,arg);
        }
        else if ( map_string.compare(0,1,"#") != 0 ) //no warings for commented-out input
            printf("\nWARNING: INPUT SCRIPT PARSING: instruction '%s' is not meaningful, hence corresponding class cannot be found. Will skip this instruction. \n \n",
                   map_string.c_str());
    }

  private:

    template<typename U>
    U* add(U *ptr,std::string map_string)
    {
      list_.push_back(static_cast<c3poBaseInterface*>(ptr));
      map_[map_string] = list_.back();
      if(ptr != list_.back())
        printf("ASSERTION FAILED\n");
      return static_cast<U*>(list_.back());
      //return ptr;
    }

    std::vector<c3poBaseInterface*> list_;
    mutable std::map<std::string,c3poBaseInterface*> map_;
};

} //end c3po_NS

#endif

