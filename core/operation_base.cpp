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


#include "operation_base.h"
#include "string.h"
#include "comm.h"
#include "output.h"
#include "error.h"
#include "input.h"
#include <fstream>

using namespace C3PO_NS;
using std::ifstream;

/* ----------------------------------------------------------------------
   ModelBase Constructor
------------------------------------------------------------------------- */

OperationBase::OperationBase(c3po *ptr,const char *name) : 
    c3poBaseAccessible(ptr),
    name_(0),
    scope_(0),
    operationID_(-1)  //id must be set by specific operation!
{
    int n = strlen(name) + 1;
    name_ = new char[n];
    strcpy(name_,name);

    n = strlen(name) + 10 + 1;
    scope_ = new char[n];
    strcpy(scope_,"operation_");
    strcat(scope_,name);
   
   // char jsonfile[200];
 
   // if (comm().is_proc_0())
    //{
       if(input().verbose()) 
        printf("\nOperation with name %s initialized. \n", name_);
       /* sprintf(jsonfile,"/%s.json",scope_);

        std::ifstream arch(jsonfile, std::ios::in);
        if (!arch)
            error().throw_error_one(FLERR,"can not open config file ",jsonfile);

        output().write_screen_one("Parsing operation config file");
        output().write_log_one("Parsing operation config file");

    //}
    
   
     
    if(comm().is_parallel())
        error().throw_error_all(FLERR,"TODO: need to communicate settings in parallel in Input:: function");
*/
}

OperationBase::~OperationBase()
{
    delete name_;
}

