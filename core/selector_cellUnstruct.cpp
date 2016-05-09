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


#include "selector_cellUnstruct.h"
#include "selector_container.h"
#include "comm.h"
#include "mesh.h"
#include "output.h"
#include <cmath>
#include <algorithm>
//Performance test


#include "timer.h"


using namespace C3PO_NS;

/* ----------------------------------------------------------------------
   SelectorCellUnstruct Constructor
------------------------------------------------------------------------- */

SelectorCellUnstruct::SelectorCellUnstruct(c3po *ptr,const char *name)
:
SelectorBase(ptr,name)
{
  cellList_=new double[3*comm().nprocs()];
  max_=new double[3*comm().nprocs()];
  min_=new double[3*comm().nprocs()];
  tolerance_=c3po_ptr()->meshFilterWidthTolerance();
 
}


SelectorCellUnstruct::~SelectorCellUnstruct()
{
  delete cellList_;
  delete max_;
  delete min_;
 
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int SelectorCellUnstruct::findNearestCell(double x, double y, double z)
{
 //Check if the point is inside this sub-domain
 if( ! checkPosition(x, mesh().dLocalMax()[0], mesh().dLocalMin()[0] ) ) return -1;
 if( ! checkPosition(y, mesh().dLocalMax()[1], mesh().dLocalMin()[1] ) ) return -1;
 if( ! checkPosition(z, mesh().dLocalMax()[2], mesh().dLocalMin()[2] ) ) return -1;

 //Start a linear search algorithm
 
 int cellId_=0;
 
 int NofCells_=mesh().NofCells();
 
 double distance_ =  (  mesh().CellCentre(0)[0] - x) * (  mesh().CellCentre(0)[0] - x)
                    + (  mesh().CellCentre(0)[1] - y) * (  mesh().CellCentre(0)[1] - y)
                    + (  mesh().CellCentre(0)[2] - z) * (  mesh().CellCentre(0)[2] - z);
 
 for(int cell=1;cell<NofCells_;cell++)
 {
  
  double distance_new_ =  ( mesh().CellCentre(cell)[0] - x) * ( mesh().CellCentre(cell)[0] - x)
                        + ( mesh().CellCentre(cell)[1] - y) * ( mesh().CellCentre(cell)[1] - y)
                        + ( mesh().CellCentre(cell)[2] - z) * ( mesh().CellCentre(cell)[2] - z);
  
  if( distance_new_ < distance_)
  {
   cellId_=cell;
   distance_=distance_new_;
  }
 
 }
 
 return cellId_;

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void SelectorCellUnstruct::begin_of_step()
{
  selectorContainer().resetCellsInFilter();
  RunSelector();
  
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void SelectorCellUnstruct::RunSelector()
{

 
 fillArrays(); 
 boundaryCorrections();
 
 timer().stamp();
 
 
 int nprocs=comm().nprocs();
// double temp_[nprocs];
//   double vol[nprocs];
   

 
 double r_ = (*selectorContainer().getFilterSize(0)) * (*selectorContainer().getFilterSize(0));
 double delta_[3];
 for(int j=0;j<3;j++)
  delta_[j] = mesh().dGlobalMax()[j] - mesh().dGlobalMin()[j];
 
 
 //Create the list 
 
 int NofCells_=mesh().NofCells();
 
 for(int p=0; p<nprocs; p++)
 {
  
  int count_=0;
//  vol[p]=0.0;
  bool loop=true;
  bool loop2=false;
  int index_=3*p;
  
 //Check if in "dummy mode"

  if((currentCell_>= mesh().MaxNofCellsProc()[p])  )
  {
   loop=false;
  }
  
  if(selectorContainer().particleBased() && (currentPar_>= dataStorage().getNofParticlesProc(p) ))
  {
   loop=false;
  }
  
   if(!inside_[p]) loop = false;
  
 //Check if within the boundaries. This can save computational time but pose some constraints on the domain decomposition technique.
 //check for periodic and additional check for boundaries
 if(loop)
 {
  for(int j=0;j<3;j++)
   if(max_[index_+j] > min_[index_+j] )
   {
    periodic_[j]=false;
    if( ( (mesh().dLocalMin()[j] ) > max_[index_ +j] ) || ( ( mesh().dLocalMax()[j]  ) < min_[index_ +j] )) 
    {
     loop=false;
     break;
    }
   
    if(j==0)      checkX=&SelectorCellUnstruct::checkPosition;
    else if(j==1) checkY=&SelectorCellUnstruct::checkPosition;
    else if(j==2) checkZ=&SelectorCellUnstruct::checkPosition;
   }
   else 
   {
    periodic_[j]=true;
    if( ( ( mesh().dLocalMin()[j] )   >  max_[index_ +j] ) && ( ( mesh().dLocalMax()[j] ) < min_[index_ +j] ) ) 
    {
     loop=false;
     break;
    }
    if(j==0)      checkX=&SelectorCellUnstruct::checkPositionPeriodic; 
    else if(j==1) checkY=&SelectorCellUnstruct::checkPositionPeriodic;
    else if(j==2) checkZ=&SelectorCellUnstruct::checkPositionPeriodic;
   } 
 
   if(selectorContainer().filterType()==0) checkR=&SelectorCellUnstruct::checkPositionSphereDummy;
   else if(selectorContainer().filterType()==1) checkR=&SelectorCellUnstruct::checkPositionPeriodicSphere;
  
  //Check if filtersize is bigger than the domain (that would save computational time...but just works for cartesian filters)
 
  if(selectorContainer().filterType()==0)
  {
   bool check[3];
   for(int i=0;i<3;i++)
   {
    if(
       (
        !periodic_[i]                                 &
        (mesh().dLocalMin()[i]  > min_[index_ + i] ) &&
        (mesh().dLocalMax()[i]  < max_[index_ + i] ) &&
        loop
       )
      
       ||
       
       (
         periodic_[i]                                &&
        (mesh().dLocalMax()[i]  < max_[index_ + i] ) &&
        loop
       )
      
       ||
     
       (
         periodic_[i]                                &&
        (mesh().dLocalMin()[i]  > min_[index_ + i] ) &&
        loop     
       )
    
     ) check[i]=true;
    else check[i]=false;
 
   }
 
   if( check[0] && check[1] && check[2])
   {
    loop=false;
    loop2=true;
   }
  } 

  if(input().verbose())
  {

    std::cout << "\n Processor: " << comm().me() << " select for proc: " << p
              << " loop: " << loop << " loop2: " <<loop2;
  }
 }
 //Now the real loop starts
 if(loop) 
  for(int cell=0; cell<NofCells_; cell++)
  {
   
   
   if( ! (
          (this->*checkX)( 
                        mesh().CellCentre(cell)[0] ,
                          max_[index_],
                          min_[index_]
                      ) 
          ) )     continue;
  
   if( ! (
          (this->*checkY)( 
                          mesh().CellCentre(cell)[1] ,
                          max_[index_+1],
                          min_[index_+1]
                      ) 
          ) )     continue;
   
   if( ! (
          (this->*checkZ)( 
                          mesh().CellCentre(cell)[2] ,
                          max_[index_+2],
                          min_[index_+2]
                      ) 
          ) )     continue;
   
   if( ! (
          (this->*checkR)(
                           mesh().CellCentre(cell)[0],
                           mesh().CellCentre(cell)[1],
                           mesh().CellCentre(cell)[2],
                          r_,
                          index_,
                          &delta_[0]
                         )
           ) )  continue;
   
   count_++;
   selectorContainer().addCellInFilter(cell); 
   //vol[p] += *(mesh().CellVol(cell));
   
  }
 else if(loop2) //Every cell is added to filter!
  for(int cell=0; cell<NofCells_; cell++)
  {
   
   count_++;
   selectorContainer().addCellInFilter(cell); 
  // vol[p] += *(mesh().CellVol(cell));

   
  }
 
  if(input().verbose())
 {
   
  //  std::cout << "\n Processor: " << comm().me() << " over proc: " << p <<" filterVolume: " << vol[p];
 }

 
  selectorContainer().writeKey(p,count_);
 }
 

 timer().stamp(TIME_SELECTOR);
 timer().stamp();
 //  MPI_Allreduce(&vol, &temp_, nprocs, MPI_DOUBLE,
 //              MPI_SUM, MPI_COMM_WORLD);
  
 timer().stamp(TIME_SELECTOR_MPI);
 
  //selectorContainer(). addFilterVolume(vol[comm().me()]);
 

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void SelectorCellUnstruct::fillArrays()
{
 
 timer().stamp();

 if(selectorContainer().particleBased())
 {  

  currentPar_=selectorContainer().currentParticle();
  
  double currParCoord_[3];
  int nprocs=comm().nprocs();
  
  //Every processor has to know the coordinates of every current particle
  if(currentPar_ <  dataStorage().getNofParticlesProc(comm().me()))
  {
   for(int i=0; i<3;i++)
    currParCoord_[i]= dataStorage().getParticle(currentPar_)->getpos()[i]; 
  }
  else
  { 
   for(int i=0; i<3;i++)
    currParCoord_[i]=  (mesh().dMin()[comm().me()*3 +i] + mesh().dMax()[comm().me()*3 + i])/2;
  }
  
 
  timer().stamp(TIME_SELECTOR);
  timer().stamp();
  
//  MPI_Barrier(MPI_COMM_WORLD);  
  MPI_Allgather(currParCoord_, 3, MPI_DOUBLE,  cellList_, 3, MPI_DOUBLE,MPI_COMM_WORLD);
  
  timer().stamp(TIME_SELECTOR_MPI);
  timer().stamp();
  
  //The absolute box size for every cell is calculated
  for(int p=0;p<nprocs;p++)
  {
  
   inside_[p]=true;
  
   //Check when the filter is selective
   if(selectorContainer().selectiveFilter())
   {
   
     for(int i=0;i<3;i++)
     {
      if(  cellList_[3*p+i]  > selectorContainer().maxToFilter()[i]
         || 
          cellList_[3*p+i] < selectorContainer().minToFilter()[i]
                
       )
      inside_[p]=false;
    
     }
    }
  
   
   if(currentPar_ < dataStorage().getNofParticlesProc(p) && inside_[p])
   {
    for(int j=0;j<3;j++)
    {
     max_[3*p+j] = cellList_[3*p+j] + *(selectorContainer().getFilterSize(j));
     min_[3*p+j] = cellList_[3*p+j] - *(selectorContainer().getFilterSize(j));   
    }
   }
   else
   {
    inside_[p]=false;
   }
  }
  timer().stamp(TIME_SELECTOR);
  
 }
 else
 {
  
  currentCell_= *(selectorContainer().currentCell());

  double currCellCoord_[3];
  int nprocs=comm().nprocs();
  
  //Every processor has to know the coordinates of every current cell
  if(currentCell_ < mesh().MaxNofCellsProc()[comm().me()])
  {
   for(int i=0; i<3;i++)
    currCellCoord_[i]= mesh().CellCentre(currentCell_)[i];;
  }
  else
  { 
   for(int i=0; i<3;i++)
    currCellCoord_[i]=  mesh().dMin()[comm().me()*3 +i] + mesh().dMax()[comm().me()*3 + i]/2;    
  }
  
  
  timer().stamp(TIME_SELECTOR);
  timer().stamp();
  
//  MPI_Barrier(MPI_COMM_WORLD);  
  MPI_Allgather(currCellCoord_, 3, MPI_DOUBLE,  cellList_, 3, MPI_DOUBLE,MPI_COMM_WORLD);
  
  timer().stamp(TIME_SELECTOR_MPI);
  timer().stamp();
  
  //The absolute box size for every cell is calculated
  for(int p=0;p<nprocs;p++)
  {
   
   inside_[p]=true;
   
   //Check when the filter is selective
   if(selectorContainer().selectiveFilter())
   {
   
     for(int i=0;i<3;i++)
     {
      if(  cellList_[3*p+i]  > selectorContainer().maxToFilter()[i]
         || 
          cellList_[3*p+i] < selectorContainer().minToFilter()[i]
                
       )
      inside_[p]=false;
    
     }
    }
   
   
   
   if(currentCell_ < mesh().MaxNofCellsProc()[p] && inside_[p])
   {
    for(int j=0;j<3;j++)
    {
     max_[3*p+j] = cellList_[3*p+j] + *(selectorContainer().getFilterSize(j));
     min_[3*p+j] = cellList_[3*p+j] - *(selectorContainer().getFilterSize(j));   
    }
   }
   else
   {
    for(int j=0;j<3;j++)
    {
     max_[3*p+j] = cellList_[3*p+j];
     min_[3*p+j] = cellList_[3*p+j];   
    }  
   }
  }
  timer().stamp(TIME_SELECTOR);
  
  
 }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void SelectorCellUnstruct::boundaryCorrections()
{
 timer().stamp();
 
  int nprocs=comm().nprocs();
  for(int p=0; p<nprocs;p++)
  {
   if(!inside_[p]) continue;
   for(int j=0;j<3;j++)
   {
    //If wall, go ahead
    if(input().getBound()[j]==-1) continue;
   
    double delta_ =  mesh().dGlobalMax()[j] - mesh().dGlobalMin()[j];
    
    if ( (max_[3*p+j] - min_[3*p+j]) > delta_) //Do not filter more than the entire domain
    {
     max_[3*p+j]=mesh().dGlobalMax()[j];
     min_[3*p+j]=mesh().dGlobalMin()[j];
    }
    else // adjust values to account for periodicity
    {
     if (max_[3*p+j] > mesh().dGlobalMax()[j] )  max_[3*p+j] -= delta_ ;
     if (min_[3*p+j] < mesh().dGlobalMin()[j] )  min_[3*p+j] += delta_ ;
    } 
   }
  }
timer().stamp(TIME_SELECTOR);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
bool SelectorCellUnstruct::checkPosition(double posA, double max, double min )
{
 
 if((posA + tolerance_) < max && (posA - tolerance_) > min) return true;
 else return false;
 
}
bool SelectorCellUnstruct::checkPositionPeriodic(double posA, double max, double min )
{
 if((posA - tolerance_) < max || (posA + tolerance_) > min) return true;
 else return false;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
bool SelectorCellUnstruct::checkPositionSphere(double i, double j, double k, double r_, int index_, double * delta_)
{
 
    double distance_ =   (i-cellList_[index_] )*(i-cellList_[index_] )
                       + (j-cellList_[index_+1] )*(j-cellList_[index_+1] )
                       + (k-cellList_[index_+2] )*(k-cellList_[index_+2] );
            
     if( distance_< r_) return true; 
      
 
 return false;

}


bool SelectorCellUnstruct::checkPositionPeriodicSphere(double i, double j, double k, double r_ ,int index_, double * delta_)
{
 
 if(checkPositionSphere(i,j,k, r_ , index_, delta_)) return true;
 
 double distancex_;
 double distancey_;
 double distancez_;
 
          
 for(int signx_=-1;signx_<2;signx_++)
 {
  if(periodic_[0]) distancex_= (i-cellList_[index_] + signx_ * delta_[0]  )*(i-cellList_[index_] + signx_ * delta_[0]);
  else
  {
   signx_=1;
   distancex_= (i-cellList_[index_]  )*(i-cellList_[index_] );
  }
  for(int signy_=-1;signy_<2;signy_++)
  {
   if(periodic_[1]) distancey_= (j-cellList_[index_+1]+ signy_ * delta_[1]  )*(j-cellList_[index_+1] + signy_ * delta_[1]);
   else
   {
    signy_=1;
    distancey_= (j-cellList_[index_+1] )*(j-cellList_[index_+1] );
   }
   for(int signz_=-1;signz_<2;signz_++)
   {
    if(periodic_[2]) distancez_= (k-cellList_[index_+2] + signz_ * delta_[2]  )*(k-cellList_[index_+2] + signz_ * delta_[2]);
    else
    {
     distancez_= (k-cellList_[index_+2])*(k-cellList_[index_+2]);
     signz_=1;
    }
    if( distancex_ + distancey_ + distancez_ < r_) return true; 
    
   }
  }  
 }
 
 return false;

}


