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
# include "operation_filtering.h"
#include "operation_container.h"
#include "selector_container.h"
#include "filter_base.h"
#include "string.h"
#include "mesh.h"
#include "comm.h"
#include "output.h"
#include "error.h"
#include <fstream>
#include <cmath>
#include <string>

#include "memory_ns.h"

#include "timer.h"

#define TOLERANCE_ALPHAVAL 1e-10

using namespace C3PO_NS;


/* * * * * * * * * * * * * * * * * * * Divergent Algorithm for data filtering * * * * * * * * * * * * * * * * * * * * */ 

void OperationFiltering::divergentAlgorithm()
{
  timer().stamp();

  //define some stuff
  int    selCell=*(selectorContainer().currentCell());
  int    nprocs_=comm().nprocs();
  int*   key=selectorContainer().getKey();
  double alphaValue;
  double* selCellCenters = selectorContainer().currentCellCenters();
  
  //tmp fields needed for variance calculations  
  double meanOldStorage[SF_to_filter];
  double meanNewStorage[SF_to_filter];
  double meanOldStorageVec[VF_to_filter][3];
  double meanNewStorageVec[VF_to_filter][3];

  std::vector<int>* cellList=selectorContainer().getCellsInFilter();
  
   bool haveCell = selectorContainer().haveCell();
   
   double tmpFields_[totFields_];
 
  
   if(!haveCell) selCell=-1;
  
  //Fill temporary containers with data
   if(haveCell) tmpFields_[totFields_-1] = calculateWeights( selCell);  
   else tmpFields_[totFields_-1]=0.0;

  for (int i=0;i<SF_to_filter;i++)
  {   
    if(haveCell)     tmpFields_[i] = dataStorage().SF(RMASF_[i],selCell);  
    else                tmpFields_[i] = 0.0;

  }
  
  for (int i=0;i<VF_to_filter;i++)
  {   
    if(haveCell) for (int coord=0;coord<3;coord++)       tmpFields_[SF_to_filter+3*i+coord] = dataStorage().VF(RMAVF_[i],coord,selCell);
    else            for (int coord=0;coord<3;coord++)    tmpFields_[SF_to_filter+3*i+coord] = 0.0;

  }
  
  timer().stamp(TIME_FILTER_MPI);
  MPI_Allgather(&tmpFields_[0],totFields_, MPI_DOUBLE,&fieldsPerProc_[0][0], totFields_, MPI_DOUBLE,MPI_COMM_WORLD);
  timer().stamp();

  //Perform main filtering loop
  int point=0;
  for (int k=0;k<nprocs_;k++)
  { 
        for (int j=point;j<key[k]+point;j++) 
        {
            
            
            alphaValue = fieldsPerProc_[k][totFields_-1]*kernelFunction( mesh().CellCentre((*cellList)[j]),&selCellCenters[k*3]); 
            filterWeights_[(*cellList)[j]] += alphaValue;
            double Wn = filterWeights_[(*cellList)[j]] ;
           double alphaDivW = alphaValue
                                / (
                                      Wn 
                                    + *(mesh().CellVol(selCell)) * TOLERANCE_ALPHAVAL 
                                  );

                                      
            //Loop scalar fields to filter and do variance calculation
            for (int i=0;i<SF_to_filter;i++)
            { 
                 meanOldStorage[i] = dataStorage().fSF(filterSF_[i])->value()[(*cellList)[j]];
                 meanNewStorage[i] = meanOldStorage[i] 
                                   + alphaDivW 
                                   * (fieldsPerProc_[k][i] - meanOldStorage[i]);
                         
                dataStorage().fSF(filterSF_[i])->value()[(*cellList)[j]] = meanNewStorage[i];
            }

            //Loop vector fields
            for (int i=0;i<VF_to_filter;i++)
            { 
               for (int coord=0;coord<3;coord++)
               {
                 meanOldStorageVec[i][coord] = *dataStorage().fVF(filterVF_[i])->value(coord,(*cellList)[j]);
                 meanNewStorageVec[i][coord] = meanOldStorageVec[i][coord]
                         + alphaDivW 
                         * (fieldsPerProc_[k][SF_to_filter+i*3+coord] - meanOldStorageVec[i][coord]);
                         
                 *dataStorage().fVF(filterVF_[i])->value(coord,(*cellList)[j]) = meanNewStorageVec[i][coord];
               }
            }
            
            //Loop Variance fields:
            if(!computeVariance_)            
                continue;


            double varianceOld;
            double varianceNew;

            //Loop SCALAR fields and do variance calculation
            for(int kVar=0;kVar<SF_for_varianceCalc_;kVar++)
            {
                 int i = SF_to_filter + kVar;
                 varianceOld = dataStorage().fSF(filterSF_[i])->value()[(*cellList)[j]]; //no separate storage, uses filter storage!

                 if(varianceHasCrossTerm_)
                 {
                   varianceNew = varianceOld 
                               + alphaValue
                               * (
                                    fieldsPerProc_[k][filterSFVarianceValueID_[kVar]]       //xi,n
                                   *fieldsPerProc_[k][filterSFVarianceValueSecondID_[kVar]]  //yi,n
                                  - meanOldStorage[filterSFVarianceValueID_[kVar]]                  //\overbar{xi,n-1]
                                   *meanOldStorage[filterSFVarianceValueSecondID_[kVar]]            //\overbar{yi,n-1]
                                 )
                               + Wn
                               *(
                                    meanOldStorage[filterSFVarianceValueID_[kVar]]                  //\overbar{xi,n-1]
                                   *meanOldStorage[filterSFVarianceValueSecondID_[kVar]]            //\overbar{yi,n-1]
                                  - meanNewStorage[filterSFVarianceValueID_[kVar]]                  //\overbar{xi,n]
                                   *meanNewStorage[filterSFVarianceValueSecondID_[kVar]]            //\overbar{yi,n]
                                );
                 }
                 else
                 {
                   varianceNew = varianceOld 
                               + alphaValue
                               * (
                                  fieldsPerProc_[k][filterSFVarianceValueID_[kVar]]
                                - meanOldStorage[filterSFVarianceValueID_[kVar]]
                                 )
                               * (
                                   fieldsPerProc_[k][filterSFVarianceValueID_[kVar]] 
                                - meanNewStorage[filterSFVarianceValueID_[kVar]]
                                 );
                 }       
                dataStorage().fSF(filterSF_[i])->value()[(*cellList)[j]] = varianceNew;
            }


            //Loop VECTOR fields and do variance calculation
            for(int kVar=0;kVar<VF_for_varianceCalc_;kVar++)
            {
               int firstCoord[3]  = {0,1,2};
               int secondCoord[3] = {0,1,2};
               if(filterVFVarianceComputeOffDiagonal_[kVar])
               {
                    firstCoord [1]=0;firstCoord [2]=1;
                    secondCoord[0]=1;secondCoord[1]=2;
               }
               int i = VF_to_filter + kVar;
               for (int coord=0;coord<3;coord++)
               {

                 varianceOld  = *dataStorage().fVF(filterVF_[i])->value(coord,(*cellList)[j]);
    
                 if(evaluateVarianceVectorScalarMixed_[kVar])
                 {
                   varianceNew = varianceOld 
                               + alphaValue
                               * (
                                    fieldsPerProc_[k][SF_to_filter+filterVFVarianceValueID_[kVar]*3 + coord]          //xi,n
                                   *fieldsPerProc_[k][filterSFVarianceValueVectorScalarMixed_[kVar]]                  //yi,n
                                  - meanOldStorageVec[filterVFVarianceValueID_[kVar]][coord]                          //\overbar{xi,n-1]
                                   *meanOldStorageVec[filterVFVarianceValueSecondID_[kVar]][coord]                    //\overbar{yi,n-1]
                                 )
                               + Wn
                               * (
                                    meanOldStorageVec[filterVFVarianceValueID_[kVar]][coord]                    //\overbar{xi,n-1]
                                   *meanOldStorageVec[filterVFVarianceValueSecondID_[kVar]][coord]              //\overbar{yi,n-1]
                                  - meanNewStorageVec[filterVFVarianceValueID_[kVar]][coord]                    //\overbar{xi,n]
                                   *meanNewStorageVec[filterVFVarianceValueSecondID_[kVar]][coord]              //\overbar{yi,n]
                                 );
                 }
                 else if(varianceHasCrossTerm_ || filterVFVarianceComputeOffDiagonal_[kVar] )
                 {
                   varianceNew = varianceOld 
                               + alphaValue
                               * (
                                    fieldsPerProc_[k][SF_to_filter+filterVFVarianceValueID_[kVar]*3 + firstCoord[coord]]         //xi,n
                                   *fieldsPerProc_[k][SF_to_filter+filterVFVarianceValueSecondID_[kVar]*3 + secondCoord[coord]]   //yi,n
                                  - meanOldStorageVec[filterVFVarianceValueID_[kVar]][coord]                         //\overbar{xi,n-1]
                                   *meanOldStorageVec[filterVFVarianceValueSecondID_[kVar]][coord]                   //\overbar{yi,n-1]
                                 )
                               + Wn
                               * (
                                    meanOldStorageVec[filterVFVarianceValueID_[kVar]][coord]                   //\overbar{xi,n-1]
                                   *meanOldStorageVec[filterVFVarianceValueSecondID_[kVar]][coord]             //\overbar{yi,n-1]
                                  - meanNewStorageVec[filterVFVarianceValueID_[kVar]][coord]                   //\overbar{xi,n]
                                   *meanNewStorageVec[filterVFVarianceValueSecondID_[kVar]][coord]             //\overbar{yi,n]
                 
                                 );
                 }
                 else
                 {
                   varianceNew = varianceOld 
                               + alphaValue
                               * (
                                  fieldsPerProc_[k][SF_to_filter+filterVFVarianceValueID_[kVar]*3 + coord]
                                - meanOldStorageVec[filterVFVarianceValueID_[kVar]][coord]
                                 )
                               * (
                                  fieldsPerProc_[k][SF_to_filter+filterVFVarianceValueID_[kVar]*3 + coord] 
                                - meanNewStorageVec[filterVFVarianceValueID_[kVar]][coord]
                                 );
                  }       
                
                 *dataStorage().fVF(filterVF_[i])->value(coord,(*cellList)[j]) = varianceNew;

//                  if(coord==0)
//                  {
//                     std::cout << "varianceNew: " << varianceNew << ", varianceOld: " << varianceOld 
//                               << ", alphaValue: " << alphaValue << ", Wn: " << Wn 
//                               << ", value1: " << valueContainerVectors_[filterVFVarianceValueID_[kVar]][k*3+firstCoord[coord]] 
//                               << ", value2: " << valueContainerVectors_[filterVFVarianceValueSecondID_[kVar]][k*3+secondCoord[coord]]  
//                               << ", meanOldStorageVec[" << filterVFVarianceValueID_[kVar] << "]: " << meanOldStorageVec[filterVFVarianceValueID_[kVar]][coord]
//                               << ", meanNewStorageVec[" << filterVFVarianceValueID_[kVar] << "]: " << meanNewStorageVec[filterVFVarianceValueID_[kVar]][coord]  
//                               << " \n"; 
//                  }
               }
            }
            
        }
        point+=key[k];

  }

  timer().stamp(TIME_FILTER);
// END MAIN FILTER LOOP
          
}
/*---------------------------------------------------------------------*/
void OperationFiltering::endDivergentAlgorithm()
{

  for (int cell=0;cell<mesh().NofCells();cell++) 
  {
       double alphaValue=filterWeights_[cell];

     if(!computeVariance_)
        continue;

     //Finalize variance calculation - SCALARS
     for(int kVar=0;kVar<SF_for_varianceCalc_;kVar++)
     {
          int i = SF_to_filter + kVar;
          if(std::abs(alphaValue)>(TOLERANCE_ALPHAVAL*(*mesh().CellVol(cell)))) //must multply with cell volume!
               dataStorage().fSF(filterSF_[i])->value()[cell] /=alphaValue;
          else
               dataStorage().fSF(filterSF_[i])->value()[cell] = 0.0;
     }

     //Finalize variance calculation - VECTORS
     for(int kVar=0;kVar<VF_for_varianceCalc_;kVar++)
     {
          int i = VF_to_filter + kVar;
          if(std::fabs(alphaValue)>(TOLERANCE_ALPHAVAL*(*mesh().CellVol(cell))))  //must multply with cell volume!
            for (int coord=0;coord<3;coord++)
               *dataStorage().fVF(filterVF_[i])->value(coord,cell) /= alphaValue;
          else
           for (int coord=0;coord<3;coord++)
               *dataStorage().fVF(filterVF_[i])->value(coord,cell) = 0.0;
     }
  }

}
/*------------------------------------------------------------------------*/
