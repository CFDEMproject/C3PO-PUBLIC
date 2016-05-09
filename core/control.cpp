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


#include "control.h"
#include "stdlib.h"
#include "input.h"
#include "output.h"

using namespace C3PO_NS;

/* ----------------------------------------------------------------------
   Control Constructor
------------------------------------------------------------------------- */

Control::Control(c3po *ptr) 
: 
    c3poBase(ptr),
    currStep_(0)
{
}

Control::~Control()
{
}

/* ----------------------------------------------------------------------
    Initialization
------------------------------------------------------------------------- */

void Control::run(int nStep)
{
    

}

// *************************************************************
void Control::parse_command(int narg,char const* const* arg)
{
 
  int endTime(0.0);
 // int timeStep(1.0);
    

  // parse optional args
  int iarg = 0;

/*  if (begin_pascal_init)
  { 
        pascal_ptr()->init();
     	printf("\nPaScal initialized!\n\n");
     	begin_pascal_init = false;
  }
  else 
  {

  }
*/
  while (iarg < narg) 
  {
       if (strcmp(arg[iarg],"run")==0) 
       {
            endTime=atoi(arg[iarg+1]);

            //Run the simulation
            output().write_screen_one("Control is firing up C3PO.");
            run(endTime);
            output().write_screen_one("C3PO run done.");
            iarg++;
      }
      else  iarg++;
  };

}

