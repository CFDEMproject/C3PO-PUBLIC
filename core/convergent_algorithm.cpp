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

/* * * * * * * * * * * * * * * * * * * Convergent Algorithm for data filtering * * * * * * * * * * * * * * * * * * * * */ 
void OperationFiltering::convergentAlgorithm()
{

 
 timer().stamp();

  //define some stuff
  int selPar=selectorContainer().currentParticle();
  int nprocs_=comm().nprocs();
  std::vector<int>* cellList=selectorContainer().getCellsInFilter();
  int* key = selectorContainer().getKey();
  int totvec= 3*(VF_to_filter);
  double localAlpha_;
  int point=0;
  double* selCellCenters = selectorContainer().currentCellCenters();
  
  bool haveParticle=selectorContainer().haveParticle();
  
  if(!haveParticle) selPar=-1;
  else
  {
    //dataStorage().getParticle(selPar)->setFilterVolume( *selectorContainer().filterVolume() );
  }
  setVectorToZero(&fieldsPerProc_[0][0],totFields_*nprocs_);

  //Calculating partial sum
  for (int k=0;k<nprocs_;k++)
  { 
   
      
    for (int s=point;s<key[k]+point;s++) 
    {
     
     localAlpha_= calculateWeights( (*cellList)[s]) * kernelFunction( mesh().CellCentre((*cellList)[s]),&selCellCenters[k*3]);
       
     for(int i=0;i<VF_to_filter;i++)
      for (int coord=0;coord<3;coord++)
       fieldsPerProc_[k][3*i+coord] += localAlpha_ * dataStorage().VF(RMAVF_[i],coord,(*cellList)[s]);
      
    
     for(int i=0;i<SF_to_filter;i++)  
      fieldsPerProc_[k][i+totvec] += localAlpha_ * dataStorage().SF(RMASF_[i],(*cellList)[s]);
     
     fieldsPerProc_[k][totFields_-1]+= localAlpha_;
     
    }
    

   point+=key[k];

  }
  
  timer().stamp(TIME_FILTER);
  timer().stamp();
  
   MPI_Allreduce(&fieldsPerProc_[0][0],&tmpData_[0][0], totFields_*nprocs_, MPI_DOUBLE, MPI_SUM ,MPI_COMM_WORLD); 
   
  
    
  timer().stamp(TIME_FILTER_MPI);
  timer().stamp();
  
 if(haveParticle && tmpData_[comm().me()][totFields_-1]>TOLERANCE_ALPHAVAL*(*selectorContainer().filterVolume())) 
 {
 //Finish Calculation

 //Send the filter volume to particle
   dataStorage().getParticle(selPar)->setFilterVolume( tmpData_[comm().me()][totFields_-1] ); 

   for(int i=0;i<VF_to_filter;i++)
      for (int coord=0;coord<3;coord++)
        dataStorage().getParticle(selPar)->filteredVector(filterVF_[i])[coord] = tmpData_[comm().me()][3*i+coord]/tmpData_[comm().me()][totFields_-1];
        
  
   for(int i=0;i<SF_to_filter;i++)  
    *(dataStorage().getParticle(selPar)->filteredScalar(filterSF_[i])) =  tmpData_[comm().me()][i+totvec]/tmpData_[comm().me()][totFields_-1]; 
  }  
       

 timer().stamp(TIME_FILTER);
    
  

    
 
  if(computeVariance_)            
  { 
    point=0;
    
    setVectorToZero(&VTmp_[0][0],totVar_*nprocs_);
    int varVecs=(VF_for_varianceCalc_*3);
    for (int k=0;k<nprocs_;k++)
    {
     if(tmpData_[k][totFields_-1]>TOLERANCE_ALPHAVAL)
     {
      for (int s=point;s<key[k]+point;s++) 
      {
       localAlpha_=  calculateWeights( (*cellList)[s])*kernelFunction( mesh().CellCentre((*cellList)[s]),&selCellCenters[k*3]);
      
       //Loop SCALAR fields and do variance calculation
        for(int kVar=0;kVar<SF_for_varianceCalc_;kVar++)
        {
     
         if(varianceHasCrossTerm_)
          VTmp_[k][kVar+varVecs] += localAlpha_                              
                          *( dataStorage().SF(RMASF_[filterSFVarianceValueID_[kVar]],(*cellList)[s]) - tmpData_[k][filterSFVarianceValueID_[kVar]+totvec]/tmpData_[k][totFields_-1])
                          *( dataStorage().SF(RMASF_[filterSFVarianceValueSecondID_[kVar]],(*cellList)[s]) -  tmpData_[k][filterSFVarianceValueSecondID_[kVar]+totvec]/tmpData_[k][totFields_-1]);
         else
          VTmp_[k][kVar+varVecs] += localAlpha_                             
                          *( dataStorage().SF(RMASF_[filterSFVarianceValueID_[kVar]],(*cellList)[s]) - tmpData_[k][filterSFVarianceValueID_[kVar]+totvec]/tmpData_[k][totFields_-1])
                          *( dataStorage().SF(RMASF_[filterSFVarianceValueID_[kVar]],(*cellList)[s]) - tmpData_[k][filterSFVarianceValueID_[kVar]+totvec]/tmpData_[k][totFields_-1]);
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
  
         for (int coord=0;coord<3;coord++)
         {

          if(evaluateVarianceVectorScalarMixed_[kVar])
             VTmp_[k][3*kVar+coord] += localAlpha_                              
                                      *( dataStorage().VF(RMAVF_[filterVFVarianceValueID_[kVar]],coord,(*cellList)[s]) - tmpData_[k][3*filterVFVarianceValueID_[kVar]+coord]/tmpData_[k][totFields_-1])
                                     *( dataStorage().VF(RMAVF_[filterVFVarianceValueSecondID_[kVar]],coord,(*cellList)[s]) - tmpData_[k][3*filterVFVarianceValueSecondID_[kVar]+coord]/tmpData_[k][totFields_-1]);
          else if(varianceHasCrossTerm_ || filterVFVarianceComputeOffDiagonal_[kVar] )
             VTmp_[k][3*kVar+coord] += localAlpha_                              
                                      *( dataStorage().VF(RMAVF_[filterVFVarianceValueID_[kVar]],firstCoord[coord],(*cellList)[s])  -  tmpData_[k][3*filterVFVarianceValueID_[kVar]+coord]/tmpData_[k][totFields_-1])
                                      *( dataStorage().VF(RMAVF_[filterVFVarianceValueSecondID_[kVar]],secondCoord[coord],(*cellList)[s]) -  tmpData_[k][3*filterVFVarianceValueID_[kVar]+coord]/tmpData_[k][totFields_-1]);
          else
             VTmp_[k][3*kVar+coord] += localAlpha_                              
                                      *( dataStorage().VF(RMAVF_[filterVFVarianceValueID_[kVar]],coord,(*cellList)[s]) -  tmpData_[k][3*filterVFVarianceValueID_[kVar]+coord]/tmpData_[k][totFields_-1])
                                      *( dataStorage().VF(RMAVF_[filterVFVarianceValueID_[kVar]],coord,(*cellList)[s]) -  tmpData_[k][3*filterVFVarianceValueID_[kVar]+coord]/tmpData_[k][totFields_-1]);



         }
        }
      }
     }
    }
  
  timer().stamp(TIME_FILTER);
  timer().stamp();
  
   MPI_Allreduce(&VTmp_[0][0],&Var_[0][0], totVar_*nprocs_, MPI_DOUBLE, MPI_SUM ,MPI_COMM_WORLD); 
   
  timer().stamp(TIME_FILTER_MPI);
  timer().stamp();
  if(selPar==-1) return;
  if(tmpData_[comm().me()][totFields_-1]<TOLERANCE_ALPHAVAL) return;
  
  double finalAlpha_=tmpData_[comm().me()][totFields_-1];
 
   for(int kVar=0;kVar<VF_for_varianceCalc_;kVar++)
   { 
      int i = VF_to_filter + kVar;
      for(int coord=0;coord<3;coord++)
      {
        dataStorage().getParticle(selPar)->filteredVector(filterVF_[i])[coord] =  Var_[comm().me()][kVar+coord] /finalAlpha_;     
      }
   }
   
   for(int kVar=0;kVar<SF_for_varianceCalc_;kVar++)
   { 
      int i = SF_to_filter + kVar;
       *(dataStorage().getParticle(selPar)->filteredScalar(filterSF_[i])) =  Var_[comm().me()][(VF_for_varianceCalc_*3)+kVar] /  finalAlpha_; 
   }
  

  
  }
 
   
 
    timer().stamp(TIME_FILTER);
} 
