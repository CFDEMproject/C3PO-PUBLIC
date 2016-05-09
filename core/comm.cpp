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

#include "comm.h"
#include "input.h"
//#include "output.h"
//#include "control.h"
//#include "simulation_state.h"
#include "stdlib.h"
#include "selector_container.h"

using namespace C3PO_NS;

/* ----------------------------------------------------------------------
   Constructor / Destructor
------------------------------------------------------------------------- */

Comm::Comm(c3po *ptr,MPI_Comm &communicator) : c3poBase(ptr),
    me_(-1),
    nprocs_(-1),
    world_(communicator)
{
   // do some MPI stuff

    MPI_Comm_rank(world_,&me_);
    MPI_Comm_size(world_,&nprocs_);
   // MPI_Comm_set_errhandler(world_,MPI_ERRORS_RETURN);
  

    //MPI_Abort(world_,1);

    //printf("Comm constructor called, I am %d, Input is %d\n",this,input());

}

Comm::~Comm()
{
//    MPI_Comm_free(&world_);
}

/* ----------------------------------------------------------------------
   settings
------------------------------------------------------------------------- */

void Comm::init()
{
    //TODO check if CommModel allocated
}

/* ----------------------------------------------------------------------
   settings
------------------------------------------------------------------------- */

void Comm::parse_command(int narg,char const* const* arg)
{
//    output().write_screen_one("Comm is parsing something:");
//    output().write_screen_one(arg[0]);
}

/* ----------------------------------------------------------------------
   abort simulation, called by one proc (crash) or
   finalize simulation, called by all proc (clean termination)
------------------------------------------------------------------------- */

void Comm::abort_one() const
{
    //MPI_Abort(MPI_COMM_WORLD,1);
    MPI_Abort(world_,1);
}

void Comm::finalize_all() const
{
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    exit(1);
}



/*MPI::Win* Comm::goRMA(double* x) const
{

// MPI::Info info = MPI::Info::Create();
 // info.Set("no_locks", "true");
  
  int win_size;
 
  win_size=sizeof(x);
  
  
  MPI::Win *window= new MPI::Win;
  *window =MPI::Win::Create(x,win_size,sizeof(double),MPI_INFO_NULL, MPI_COMM_WORLD);
  
  return window; 
  
 // info.Free();
    
}

MPI::Win* Comm::goRMA(int* x) const
{

 // MPI::Info info = MPI::Info::Create();
 // info.Set("no_locks", "true");
  
  int win_size;
 
  win_size=sizeof(x);
  
  
  MPI::Win *window= new MPI::Win;
  *window =MPI::Win::Create(x,win_size,sizeof(int),MPI_INFO_NULL, MPI_COMM_WORLD);
  
  return window; 
  
 // info.Free();
    
}

*/
/* ----------------------------------------------------------------------
   settings
------------------------------------------------------------------------- */

//TODO move this to Comm model
/*
void Comm::exchange()
{
    like LAMMPS exchange

    push to buf
    delete particle
    comm buffer
    recv buffer
    add partice
}
*/
