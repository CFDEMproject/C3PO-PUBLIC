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


#include "mesh.h"
#include "input.h"
#include "style_selector.h"
#include "error.h"
#include <vector>
#include <cmath>

using namespace C3PO_NS;

/* ----------------------------------------------------------------------
   c3poMesh Constructor
------------------------------------------------------------------------- */

c3poMesh::c3poMesh(c3po *ptr) 
:
c3poBase(ptr),
meshXYZspacing_(1)
{

   CellCount_global_[0]=-1;
   CellCount_global_[1]=-1;
   CellCount_global_[2]=-1;
	CellCount_[0]=-1;
	CellCount_[1]=-1;
	CellCount_[2]=-1;
	
	cellSize_[0] =-1;

    complete_domain_max= new double[3*comm().nprocs()];
    complete_domain_min= new double[3*comm().nprocs()];
    
    complete_IJKdomain_max= new int[3*comm().nprocs()];
    complete_IJKdomain_min= new int[3*comm().nprocs()];
    
    MaxNofCellsProc_ = new int[comm().nprocs()];
}

/* ---------------------------------------------------------------------- */

c3poMesh::~c3poMesh()
{
    delete complete_domain_max;
    delete complete_domain_min;
    delete complete_IJKdomain_max;
    delete complete_IJKdomain_min;
    delete MaxNofCellsProc_;
    clearCells();

}

// ------------------------------------------------------------------------
void c3poMesh::gatherInfo() const
{
 int displ[comm().nprocs()];
 int numfrags[comm().nprocs()];
 
 if(cellSize_[0] == -1) calculateIJKfromJson();
        
 int sum = 0;
     
 for (int i = 0; i < comm().nprocs(); ++i) {
         
        displ[i] = sum;
        sum +=3; 
        numfrags[i]=3; 
        }
 
 MPI_Barrier(MPI_COMM_WORLD);  
 MPI_Allgatherv(local_domain_max, numfrags[comm().me()], MPI_DOUBLE,  complete_domain_max, numfrags, displ, MPI_DOUBLE,MPI_COMM_WORLD);
 MPI_Allgatherv(local_domain_min, numfrags[comm().me()], MPI_DOUBLE,  complete_domain_min, numfrags, displ, MPI_DOUBLE,MPI_COMM_WORLD);
 //Gather NofCells
 MPI_Allgather(&NofCells_, 1, MPI_INT,  MaxNofCellsProc_, 1, MPI_INT,MPI_COMM_WORLD);  
 
 //Find the max
 MaxNofCells_=MaxNofCellsProc_[0];
 for(int p=1;p<comm().nprocs();p++)
  if(MaxNofCells_<MaxNofCellsProc_[p]) MaxNofCells_=MaxNofCellsProc_[p];
  
 //Gather IJK info

 for(int i=0;i<3;i++)
 {
   local_ijk_max[i]= lround( ( local_domain_max[i] - global_domain_min[i])/cellSize_[i] );
   local_ijk_min[i]= lround( ( local_domain_min[i] - global_domain_min[i])/cellSize_[i] );
   CellCount_global_[i] = lround( ( global_domain_max[i] - global_domain_min[i])/cellSize_[i] );
   CellCount_[i]=local_ijk_max[i] -local_ijk_min[i];
 }
 
 MPI_Barrier(MPI_COMM_WORLD);  
 MPI_Allgatherv(local_ijk_max, numfrags[comm().me()], MPI_INT,  complete_IJKdomain_max, numfrags, displ, MPI_INT,MPI_COMM_WORLD);
 MPI_Allgatherv(local_ijk_min, numfrags[comm().me()], MPI_INT,  complete_IJKdomain_min, numfrags, displ, MPI_INT,MPI_COMM_WORLD);
 
}
// ------------------------------------------------------------------------
void c3poMesh::registerCells(double* Vvector, double* posVector, int NumberOfCells, int meshXYZspacing) 
{
 cellV_ = Vvector;
 cellCentre_ = posVector;
 NofCells_ = NumberOfCells;
 meshXYZspacing_ = meshXYZspacing;
}
// ------------------------------------------------------------------------
void c3poMesh::clearCells() const
{

/* cellV_.clear();
 for(unsigned int i=0;i<cellCentre_.size();i++)
  delete cellCentre_[i];
  
 cellCentre_.clear();
*/
}
// ------------------------------------------------------------------------
void c3poMesh::registerDomainInfo(double maxDomain[3],double minDomain[3],double maxDomainGlobal[3],double minDomainGlobal[3]) const
{
 for(int i=0; i<3;i++)
 {
  global_domain_max[i]=maxDomainGlobal[i];
  global_domain_min[i]=minDomainGlobal[i];
      
  local_domain_max[i]=maxDomain[i];
  local_domain_min[i]=minDomain[i];
 }
 gatherInfo(); 

}

// ------------------------------------------------------------------------
double c3poMesh::meshCheckTolerance() const
{
    return input().meshMainSettings()["checkTolerance"].toDouble();
}

// ------------------------------------------------------------------------
double c3poMesh::meshFilterWidthTolerance() const
{
    return input().meshMainSettings()["filterWidthTolerance"].toDouble();
}

// ------------------------------------------------------------------------
bool c3poMesh::meshVerbose() const
{
    return input().meshMainSettings()["verbose"].toBool();
}
//-------------------------------------------------------------------------
void c3poMesh::setcell(int* CellCount) const
{

 for(int i=0;i<3;i++)
  CellCount_[i] = CellCount[i];
 
}
//-------------------------------------------------------------------------
void c3poMesh::setcell_global(int* CellCount) const
{

 for(int i=0;i<3;i++)
  CellCount_global_[i] = CellCount[i];
 
}
//--------------------------------------------------------------------------
void c3poMesh::calculateIJKfromJson() const
{
 for(int i=0;i<3;i++)
  cellSize_[i] = input().cellSizefromJson(i);
}

