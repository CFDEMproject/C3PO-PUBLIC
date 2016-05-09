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


#include "selector_cellIJK.h"
#include "selector_container.h"
#include "comm.h"
#include "mesh.h"
#include "output.h"
#include <cmath>
//Performance test


#include "timer.h"


using namespace C3PO_NS;

/* ----------------------------------------------------------------------
   SelectorCellIJK Constructor
------------------------------------------------------------------------- */

SelectorCellIJK::SelectorCellIJK(c3po *ptr,const char *name)
:
SelectorBase(ptr,name)
{
  IJKlist_=new int[3*comm().nprocs()];
  cellList_=new double[3*comm().nprocs()];
  max_=new int[3*comm().nprocs()];
  min_=new int[3*comm().nprocs()];
  tolerance_=c3po_ptr()->meshFilterWidthTolerance();
  cellCoord_=new double[3*comm().nprocs()];

 
}

SelectorCellIJK::~SelectorCellIJK()
{
  delete IJKlist_;
  delete max_;
  delete min_;

 
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int SelectorCellIJK::findNearestCell(double x, double y, double z)
{
 
 int ijk_x =  lround(   ( x - mesh().dLocalMin()[0] )
                          / mesh().cellSize()[0]
                          
                    );
 
 int ijk_y =  lround(   ( y - mesh().dLocalMin()[1] )
                          / mesh().cellSize()[1]
                         
                    );
 

 int ijk_z =  lround(   ( z - mesh().dLocalMin()[2] )
                          / mesh().cellSize()[2]
                         
                    );
 
 int cellId_= ijk2CellIDOF(ijk_x, ijk_y, ijk_z); 
 
  //Check if the point is inside this sub-domain
 if (cellId_ >= mesh().NofCells() || cellId_<0) return -1;

 return cellId_;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void SelectorCellIJK::begin_of_step()
{
  selectorContainer().resetCellsInFilter();
  
  RunSelector();
  
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void SelectorCellIJK::RunSelector()
{

 fillArrays(); 
 boundaryCorrections();
 int nprocs=comm().nprocs();
 //double temp_[nprocs];
// double vol[nprocs];
 
 
 double r_ = (*selectorContainer().getFilterSize(0)) * (*selectorContainer().getFilterSize(0));
 double delta_[3];
 for(int j=0;j<3;j++)
  delta_[j] = mesh().dGlobalMax()[j] - mesh().dGlobalMin()[j];
 
 
 //Create the list 
 
 
 for(int p=0; p<nprocs; p++)
 {
  
  int count_=0;
 // vol[p]=0.0;
  bool loop=true;
  bool loop2=false;
  bool periodic[3];
  int index_=3*p;

 //Check if in "dummy mode"
 
  if((currentCell_>= mesh().MaxNofCellsProc()[p]))
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
    if( ( mesh().localIJKMin()[j]  >= max_[index_ +j] ) || (  mesh().localIJKMax()[j]  <= min_[index_ +j] )) 
    {
     loop=false;
     break;
    }
    
    if (  mesh().localIJKMin()[j] >= min_[index_ +j] )    min_[index_ +j] = mesh().localIJKMin()[j];
    if (  mesh().localIJKMax()[j] <= max_[index_ +j] )   max_[index_ +j] = mesh().localIJKMax()[j];
    periodic[j]=false;
   }
   else 
   {
    periodic[j]=true;
    
   
    if( ( ( mesh().localIJKMin()[j] )   >=  max_[index_ +j] ) && ( ( mesh().localIJKMax()[j] ) <= min_[index_ +j] ) ) 
    {
     loop=false;
     break;
    }
    
    if( mesh().localIJKMin()[j] >= min_[index_ +j] )
    {
     min_[index_ +j] = mesh().localIJKMin()[j];
     max_[index_ +j] = mesh().localIJKMax()[j];
     
     periodic[j]=false;
    } 
   
    if( mesh().localIJKMax()[j] <= max_[index_ +j] )
    {
     min_[index_ +j] = mesh().localIJKMin()[j];
     max_[index_ +j] = mesh().localIJKMax()[j];
     
     periodic[j]=false;
    } 
  
   } 
 
   if(selectorContainer().filterType()==0) checkR=&SelectorCellIJK::checkPositionSphereDummy;
   else if(selectorContainer().filterType()==1) checkR=&SelectorCellIJK::checkPositionPeriodicSphere;
  
  if(selectorContainer().filterType()==0)
  {
   bool check[3];
   for(int i=0;i<3;i++)
   {
    if(
       (
        !periodic_[i]                                 &
        (mesh().localIJKMin()[i]  > min_[index_ + i]) &&
        (mesh().localIJKMax()[i]  < max_[index_ + i]) &&
        loop
       )
      
       ||
       
       (
         periodic_[i]                                &&
        (mesh().localIJKMax()[i]  < max_[index_ + i] ) &&
        loop
       )
      
       ||
      
       (
         periodic_[i]                                &&
        (mesh().localIJKMin()[i]  > min_[index_ + i]) &&
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
 }
  
 //Entering the loop
  if(loop) 
  { 
   
   int Ncell_[2][3];
   int startId_[2][3];
  
   for(int i=0;i<2;i++) 
    for(int j=0;j<3;j++)
    {
     Ncell_[i][j]=0;
     startId_[i][j]=0;
    }
   
   
   
   //Prepare IJK info
   for(int j=0;j<3;j++)
   {
   
    if(periodic[j])
    {
     
     if( ( max_[index_ +j]   > mesh().localIJKMin()[j] )  &&  ( max_[index_ +j]  <=  mesh().localIJKMax()[j]  )   )
     {
     
      Ncell_[0][j] = max_[index_ +j] - mesh().localIJKMin()[j];
    
     }
     
     
     if( ( min_[index_ +j]   <  mesh().localIJKMax()[j] )  &&  ( min_[index_ +j]  > mesh().localIJKMin()[j] ) )
     {
      Ncell_[1][j] = mesh().localIJKMax()[j] - min_[index_ +j];
      
      startId_[1][j]= min_[index_ +j] - mesh().localIJKMin()[j];     
     }
     
    
    }
    else
    {
     Ncell_[0][j] = max_[index_ +j] - min_[index_ +j] ;
                        
     startId_[0][j]= min_[index_ +j] - mesh().localIJKMin()[j];
    } 
  
     
   }  
   
   int cell_;
  
   //This is the real IJK loop 
   for(int a=0;a<2;a++)
    for(int k=startId_[a][2];k < startId_[a][2] + Ncell_[a][2]; k++)
     for(int b=0;b<2;b++)
      for(int j=startId_[b][1];j < startId_[b][1] + Ncell_[b][1]; j++)
       for(int c=0;c<2;c++)
        for(int i=startId_[c][0];i < startId_[c][0] + Ncell_[c][0]; i++)
        {
         cell_ = ijk2CellIDOF(i,j,k);
         if(!(this->*checkR)(i,j,k,r_,index_,&delta_[0]))  continue;
        
         count_++;
         selectorContainer().addCellInFilter(cell_); 
         

        }
   }
   else if(loop2)
   {
    int NofCells_=mesh().NofCells();
    for(int cell=0; cell<NofCells_; cell++)
    {
   
     count_++;
     selectorContainer().addCellInFilter(cell); 
     
    }
   }
      
  
  //MPI_Barrier(MPI_COMM_WORLD);
  
 // vol[p] = count_ * ( *(mesh().CellVol(0)) ); 
 
  selectorContainer().writeKey(p,count_);
 }
 
 timer().stamp(TIME_SELECTOR);
 timer().stamp();
// MPI_Allreduce(&vol[0], &temp_[0], nprocs, MPI_DOUBLE,
//               MPI_SUM, MPI_COMM_WORLD);
 timer().stamp(TIME_SELECTOR_MPI);
 timer().stamp();
 // selectorContainer(). addFilterVolume(temp_[comm().me()]);

 timer().stamp(TIME_SELECTOR);

}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void SelectorCellIJK::fillArrays()
{
 timer().stamp();
 
 if(selectorContainer().particleBased())
 {  

  currentPar_=selectorContainer().currentParticle();

 // int currParCoord_[3];
  int nprocs=comm().nprocs();
  
 
 // MPI_Allgather(currParCoord_, 3, MPI_INT,  IJKlist_, 3, , MPI_INT,MPI_COMM_WORLD);
  
  double parcoordinate_[3];
  
  if(currentPar_ <  dataStorage().getNofParticlesProc(comm().me()))
  { 
   for(int i=0; i<3;i++)
    parcoordinate_[i]= dataStorage().getParticle(currentPar_)->getpos()[i]; 
  }
  else
  {
   for(int i=0; i<3;i++)
    parcoordinate_[i]= 0.0; 
  
  }
  
 // if(selectorContainer().filterType()==1) 
   MPI_Allgather(parcoordinate_, 3, MPI_DOUBLE,  cellCoord_, 3, MPI_DOUBLE,MPI_COMM_WORLD);
   
  timer().stamp(TIME_SELECTOR_MPI);
  timer().stamp();
  
 for(int p=0;p<nprocs;p++)
 {
  //Every processor has to know the coordinates of every current particle
   
   for(int i=0; i<3;i++)
    IJKlist_[p*3+i]= lround(  ( cellCoord_[p*3+i]-  mesh().dGlobalMin()[i]  ) / mesh().cellSize()[i]);
 } 

  
  //The absolute box size for every cell is calculated
  for(int p=0;p<nprocs;p++)
  {
   
   inside_[p]=true;
   
   //Check when the filter is selective
   if(selectorContainer().selectiveFilter())
   {
   
     for(int i=0;i<3;i++)
     {
      if(  IJKlist_[p*3+i]  > ( selectorContainer().maxToFilter()[i] -   mesh().dGlobalMin()[i]  )/( mesh().cellSize()[i] )
         || 
          IJKlist_[p*3+i] < ( selectorContainer().minToFilter()[i] -   mesh().dGlobalMin()[i]  )/( mesh().cellSize()[i] )
                
        )
        inside_[p]=false;
    
     }
    }
  
   
   
   if(currentPar_ <  dataStorage().getNofParticlesProc(p) && inside_[p])
   {
    for(int j=0;j<3;j++)
    {
     max_[3*p+j] = IJKlist_[3*p+j] + selectorContainer().filterWidth()[j];
     min_[3*p+j] = IJKlist_[3*p+j] - selectorContainer().filterWidth()[j] +1;   
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
  
 // int currCellIJK_[3];
  int nprocs=comm().nprocs();
 
  
  double currCellCoord_[3];
 
   //Every processor has to know the coordinates of every current cell
   if(currentCell_ < mesh().MaxNofCellsProc()[comm().me()])
   {
    for(int i=0; i<3;i++)
     currCellCoord_[i]=  mesh().CellCentre(currentCell_)[i];
   }
   else
   {
    for(int i=0; i<3;i++)
     currCellCoord_[i]=  mesh().dMin()[comm().me()*3 +i] + mesh().dMax()[comm().me()*3 + i]/2;
   }
  
  
 
 // MPI_Barrier(MPI_COMM_WORLD);  
 // MPI_Allgather(currCellIJK_, 3, MPI_INT,  IJKlist_, 3,  MPI_INT,MPI_COMM_WORLD);
  //if(selectorContainer().filterType()==1) 
   MPI_Allgather(currCellCoord_, 3, MPI_DOUBLE,  cellCoord_, 3, MPI_DOUBLE,MPI_COMM_WORLD);
 
  for(int p=0;p<nprocs;p++)
  {
  //Every processor has to know the coordinates of every current cell  
   for(int i=0; i<3;i++)
    IJKlist_[p*3+i]= lround(  (cellCoord_[p*3+i]-  mesh().dGlobalMin()[i]  ) / mesh().cellSize()[i]);
  } 

  
 
 //The absolute box size for every cell is calculated
  for(int p=0;p<nprocs;p++)
  {
   
        inside_[p]=true;
   
   //Check when the filter is selective
   if(selectorContainer().selectiveFilter())
   {
   
     for(int i=0;i<3;i++)
     {
      if(  IJKlist_[3*p+i]  > ( selectorContainer().maxToFilter()[i] -  mesh().dGlobalMin()[i] )/( mesh().cellSize()[i] )
         || 
           IJKlist_[3*p+i]  < ( selectorContainer().minToFilter()[i] -  mesh().dGlobalMin()[i] )/( mesh().cellSize()[i] )
                
       )
      inside_[p]=false;
    
     }
    }
   

  
   if(currentCell_ < mesh().MaxNofCellsProc()[p] && inside_[p])
   {
    for(int j=0;j<3;j++)
    {
     max_[3*p+j] = IJKlist_[3*p+j] + selectorContainer().filterWidth()[j];
     min_[3*p+j] = IJKlist_[3*p+j] - selectorContainer().filterWidth()[j] +1;   
    }
   }
   else
   {
    for(int j=0;j<3;j++)
    {
     max_[3*p+j] = IJKlist_[3*p+j];
     min_[3*p+j] = IJKlist_[3*p+j];   
    }  
   }
  }
 
 }
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void SelectorCellIJK::boundaryCorrections()
{
  int nprocs=comm().nprocs();
  for(int p=0; p<nprocs;p++)
  {
   if(!inside_[p]) continue;
   for(int j=0;j<3;j++)
   {
    //If wall, go ahead
    if(input().getBound()[j]==-1) continue;
   
    
    if ( (max_[3*p+j] - min_[3*p+j]) >  mesh().cell_global()[j] ) //Do not filter more than the entire domain
    {
     max_[3*p+j]= mesh().cell_global()[j];
     min_[3*p+j]= 0;
    }
    else // adjust values to account for periodicity
    {
     if (max_[3*p+j] >  mesh().cell_global()[j] )  max_[3*p+j] -=  mesh().cell_global()[j];
     if (min_[3*p+j] <  0 )                                     min_[3*p+j] +=  mesh().cell_global()[j];
    } 
   }
  }
}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
int SelectorCellIJK::ijk2CellIDOF(int i_add_correct, int j_add_correct, int k_add_correct) 
{
    int CellID_OF;          
   
    CellID_OF =   (k_add_correct)*(  mesh().cells_subdomain()[0]  * mesh().cells_subdomain()[1]   )
                + (j_add_correct)* mesh().cells_subdomain()[0] 
                +  i_add_correct;
    
    return CellID_OF;
}
/* * * * * * * * * * * * * * * * ** * * * * * * * * * * * * *  * * * *  * * */
void SelectorCellIJK::CellID2ijk(int center_Cell_ID, int* result)
{
 ////////////////////////   Cell ID to i-j-k indexing   ///////////////////////////
   
 int x=mesh().cells_subdomain()[0];  
    
 int y=mesh().cells_subdomain()[1];
 result[2] = int(center_Cell_ID / (x*y)); //estimation of z-layer
 result[1] = int((center_Cell_ID-(x*y*result[2]))/x); //estimation of y-layer
 result[0] = int((center_Cell_ID-(x*y*result[2]))-((result[1]*x))); //estimation of x-layer
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
bool SelectorCellIJK::checkPosition(double posA, double max, double min )
{
 
 if((posA + tolerance_) < max && (posA - tolerance_) > min) return true;
 else return false;
 
}
bool SelectorCellIJK::checkPositionPeriodic(double posA, double max, double min )
{
 if((posA - tolerance_) < max || (posA + tolerance_) > min) return true;
 else return false;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
bool SelectorCellIJK::checkPositionSphere(double i_ijk, double j_ijk, double k_ijk, double r_, int cell_, double * delta_)
{
 
 double i = i_ijk * mesh().cellSize()[0] + mesh().cellSize()[0]/2 + mesh().dLocalMin()[0];
 double j = j_ijk * mesh().cellSize()[1] + mesh().cellSize()[1]/2  + mesh().dLocalMin()[1];
 double k = k_ijk * mesh().cellSize()[2] + mesh().cellSize()[2]/2  + mesh().dLocalMin()[2];
 
 double distance_ =   (i-cellCoord_[cell_] )*(i-cellCoord_[cell_] )
                       + (j-cellCoord_[cell_+1] )*(j-cellCoord_[cell_+1] )
                       + (k-cellCoord_[cell_+2])*(k-cellCoord_[cell_+2] );
            
 if( distance_  - r_  < tolerance_) return true; 
      
 
 return false;

}


bool SelectorCellIJK::checkPositionPeriodicSphere(double i_ijk, double j_ijk, double k_ijk, double r_, int cell_, double * delta_)
{
 
 if(checkPositionSphere(i_ijk,j_ijk,k_ijk, r_ , cell_, delta_)) return true;
 
 double i = i_ijk * mesh().cellSize()[0] + mesh().cellSize()[0]/2 + mesh().dLocalMin()[0];
 double j = j_ijk * mesh().cellSize()[1] + mesh().cellSize()[1]/2  + mesh().dLocalMin()[1];
 double k = k_ijk * mesh().cellSize()[2] + mesh().cellSize()[2]/2  + mesh().dLocalMin()[2];
 
 double distancex_;
 double distancey_;
 double distancez_;
 
          
 for(int signx_=-1;signx_<2;signx_++)
 {
  if(periodic_[0]) distancex_= (i-cellCoord_[cell_] + signx_ * delta_[0]  )*(i-cellCoord_[cell_] + signx_ * delta_[0]);
  else
  {
   signx_=1;
   distancex_= (i-cellCoord_[cell_]  )*(i-cellCoord_[cell_] );
  }
  for(int signy_=-1;signy_<2;signy_++)
  {
   if(periodic_[1]) distancey_= (j-cellCoord_[cell_+1] + signy_ * delta_[1]  )*(j-cellCoord_[cell_+1] + signy_ * delta_[1]);
   else
   {
    signy_=1;
    distancey_= (j-cellCoord_[cell_+1] )*(j-cellCoord_[cell_+1] );
   }
   for(int signz_=-1;signz_<2;signz_++)
   {
    if(periodic_[2]) distancez_= (k-cellCoord_[cell_+2] + signz_ * (delta_[2] )  )*(k-cellCoord_[cell_+2] + signz_ * (delta_[2]));
    else
    {
     signz_=1;
     distancez_= (k-cellCoord_[cell_+2])*(k-cellCoord_[cell_+2]);
    }
    if( distancex_ + distancey_ + distancez_ + tolerance_ - r_ <tolerance_ ) return true; 
    
   }
  }  
 }
 
 return false;

}

