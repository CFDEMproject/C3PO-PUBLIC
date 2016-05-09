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

// C or Fortran style library interface to c3po
// customize by adding new functions

#include "mpi.h"
#include "stdlib.h"
#include "library.h"
#include "c3po.h"

using namespace C3PO_NS;

/* ----------------------------------------------------------------------
   create an instance of c3po and return pointer to it
   pass in command-line args and MPI communicator to run on
------------------------------------------------------------------------- */
void c3po_open(int argc, char **argv, MPI_Comm communicator, void **ptr)
{
  c3po *c3poObject = new c3po(argc,argv,communicator);
  *ptr = (void *) c3poObject;
}

/* ----------------------------------------------------------------------
   create an instance of c3po and return pointer to it
   caller doesn't know MPI communicator, so use MPI_COMM_WORLD
   intialize MPI if needed
------------------------------------------------------------------------- */

void c3po_open_no_mpi(int argc, char **argv, void **ptr)
{
  int flag;
  MPI_Initialized(&flag);

  if (!flag) {
    int argc = 0;
    char **argv = NULL;
    MPI_Init(&argc,&argv);
  }

  MPI_Comm communicator = MPI_COMM_WORLD;

  c3po * myC3po= new c3po(argc,argv,communicator);
  *ptr = (void *) myC3po;
}

