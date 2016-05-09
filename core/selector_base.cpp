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


#include "selector_base.h"
#include "string.h"
#include "comm.h"
#include "output.h"
#include "error.h"
#include <fstream>

using namespace C3PO_NS;
using std::ifstream;

/* ----------------------------------------------------------------------
   ModelBase Constructor
------------------------------------------------------------------------- */

SelectorBase::SelectorBase(c3po *ptr,const char *name) : c3poBaseAccessible(ptr),
    name_(0),
    scope_(0),
    currentCell_(-2),
    currentPar_(-2)
{
    int n = strlen(name) + 1;
    name_ = new char[n];
    strcpy(name_,name);

    n = strlen(name) + 9 + 1;
    scope_ = new char[n];
    strcpy(scope_,"selector_");
    strcat(scope_,name);
    inside_=new bool[comm().nprocs()];
 
       if(input().verbose())
        printf("\nSelector with name %s initialized. \n", name_);

}

SelectorBase::~SelectorBase()
{
    delete name_;
    delete inside_;
}

