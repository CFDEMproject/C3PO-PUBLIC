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

#include "c3po_CSV_interface.h"
#include "mpi.h"
#include <cstdio>

using namespace C3PO_NS;

int main(int argc, char **argv)
{
 
 MPI_Init(&argc,&argv);
 int nprocs_=1;
 MPI_Comm_size(MPI_COMM_WORLD,&nprocs_);
 
 c3poCSVinterface* c3po_ = new c3poCSVinterface(MPI_COMM_WORLD);
 
 c3po_->checkMesh();
 c3po_->checkParticles();
 c3po_->FLUENTrunC3PO();
 

 MPI_Barrier(MPI_COMM_WORLD);
 c3po_->clearParticles();
 
 delete c3po_;
 std::cout << "\n**************************** \n";
 std::cout << "*c3po_csv - end of program * \n";
 std::cout << "**************************** \n";
 MPI_Barrier(MPI_COMM_WORLD);
 if(nprocs_>1) MPI_Finalize();
 return 0;
}
