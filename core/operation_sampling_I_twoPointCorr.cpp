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
#include "operation_sampling_I_twoPointCorr.h"
#include "operation_container.h"
#include "data_storage.h"
#include "error.h"
#include "mesh.h"
#include "selector_container.h"
#include <fstream>
#include <sys/stat.h>
#include <iomanip>
#include <cmath>
#include "input.h"
#include "output.h"

using namespace C3PO_NS;
using namespace std;

TwoPointCorr::TwoPointCorr(c3po *ptr,const char *name) 
: 
OperationSampling(ptr,name),
component_(-1)
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
TwoPointCorr::~TwoPointCorr()
{
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /

void TwoPointCorr::begin_of_step()
{
  (this->*run)();
}

/* ----------------------------------------------------------------------*/

void TwoPointCorr::process_input(QJsonObject jsonObj)
{
 createSampleVectors(2,0,0);
        
 //Error Checks
 if(jsonObj["phaseFractionField"].isNull())
  error().throw_error_one("operation_sampling.cpp",
                                     0,
                                    "ERROR: You must specify the 'phaseFractionField' field in the Json file when using TwoPointsCorrelation \n");
 if(sampleCount_<1)
  error().throw_error_one(FLERR,"Your sampleCount is zero or not specified. \n");

 if(sampleDelta_<1e-16)
  error().throw_error_one(FLERR,"Your sampleDelta is too small or not specified. \n");
            
 if(!jsonObj["direction"].isNull())
  component_ = jsonObj["direction"].toInt();
 else
  error().throw_error_one("operation_sampling.cpp",
                                0,
                                "ERROR:You should specify a sampling direction for the two points correlation!");
         
 if(component_>2 || component_<0) 
  error().throw_error_one("operation_sampling.cpp",
                                0,
                                "ERROR:Invalid direction! Valid directions are:\n 0 (= x )\n 1 (= y )\n 2 (= z ) \n");
       
        
    
 QString Qalpha=jsonObj["phaseFractionField"].toString();
 std::string alpha=Qalpha.toUtf8().constData();
 alphaName_.assign(alpha);
       
 if((input().getSFnumber(alphaName_))==-1) dataStorage().addFieldToConvert(alphaName_);
            
 QString fN=jsonObj["fieldToSample"].toString(); 
 VFtoSample_.push_back(fN.toUtf8().constData());
       
 if(!jsonObj["filteredField"].isNull())
 {
  QString aField=jsonObj["filteredField"].toString();
  std::string aField_=aField.toUtf8().constData();
  averagedField_.assign(aField_);  
  
  if((input().getVFnumber(averagedField_))==-1) dataStorage().addFieldToConvert(averagedField_);
          
  run=&TwoPointCorr::twoPointsCorrFluct;      
 }
 else
  run=&TwoPointCorr::twoPointsCorr;  


}

//---------------------------------------------------------------------------------------------

void TwoPointCorr::twoPointsCorr()
{
  int id=dataStorage().fVFid(VFtoSample_[0]);
  if(id==-1) error().throw_error_one("OperationSampling::TwoPointsCorr()",0,"ERROR: Vector field not registered!!"); 
  
  int totcells_;
  int celli;
  int probeID;
  
  double samplePos_[3];
  
  double v1[3];
  double v2[3];
  double alpha1, alpha2;
  double sampleValue;
  
  int cellID;

  double x1_;
  
  if(lagrangian_) dataStorage().setProbes(probesName_);
  
  alpha_=dataStorage().fSFid(alphaName_);
  if(alpha_==-1) error().throw_error_one("OperationSampling::TwoPointsCorr()",0,"ERROR: phase fraction field not registered!!"); 
 
  if(lagrangian_)
   totcells_=dataStorage().numOfParticles();
  else
   totcells_=mesh().NofCells();

  int j=component_;
    
  for(int rID=0; rID<sampleCount_; rID++)
  {
      x1_= double(rID)*sampleDelta_;
     
      
      
      for (int cell=0;cell<totcells_;cell++)
      { 
       
       
       if(lagrangian_)
       {
       //calculate CellID
       
       samplePos_[0] = dataStorage().getParticle(cell)->getpos()[0];
       samplePos_[1] = dataStorage().getParticle(cell)->getpos()[1];
       samplePos_[2] = dataStorage().getParticle(cell)->getpos()[2];
       
       cellID = selectorContainer().selectCell(   samplePos_[0],
                                                  samplePos_[1],
                                                  samplePos_[2]
                                              );
       }
       else
       {
        cellID=cell; 
      
        samplePos_[0] = (mesh().CellCentre(cell)[0]);
        samplePos_[1] = (mesh().CellCentre(cell)[2]);
        samplePos_[2] = (mesh().CellCentre(cell)[1]);    
        
        if(dataStorage().SF(alpha_,cellID) > (1-1e-05) ) //no particle in there
         continue;
         
       }
               
       if(cellID>=mesh().NofCells()) error().throw_error_one("OperationSampling::TwoPointsCorr()",0,"ERROR: Problems with cellID!!"); 
      
       //Location for sampling
         
       samplePos_[j] += x1_;

       probeID= selectorContainer().selectCell(   samplePos_[0],
                                                  samplePos_[1],
                                                  samplePos_[2]
                                              );
       
       
       if(probeID>=mesh().NofCells()) 
         error().throw_error_one("OperationSampling::TwoPointsCorr()",0,"ERROR: Problems with probeID!!"); 
       
       if(probeID == -1) continue;      //TODO implement parallel TpointCorr!!!
     
    
       if(lagrangian_)
        alpha1 = *(dataStorage().getParticle(cell)->filteredScalar(alpha_));
       else
        alpha1 =dataStorage().SF(alpha_,cellID);
       
        alpha2 =dataStorage().SF(alpha_,probeID); 
    
    
        //pull out vector from storage
       
       for (int i=0;i<3;i++)
       {
        if(lagrangian_)
         v1[i]= dataStorage().getParticle(cell)->filteredVector(id)[i];
        else
         v1[i]= dataStorage().VF(id,i,cellID);
        
        v2[i]= dataStorage().VF(id,i,probeID);
       }
    
     
       // 1 - alphaAlphaUsUs 
       sampleValue = alpha1 * alpha2 * (   v1[0] * v2[0] 
                                         + v1[1] * v2[1] 
                                         + v1[2] * v2[2]
                                       );
       if( alpha2 < (1-1e-05) )
       {
          insertSample(sampleValue, x1_); //pushes the sample into the container
          if(input().verbose())
          {
            std::cout << "\n"
                 << "alpha1: "<< alpha1 << "   alpha2: " << alpha2 << "    " 
                 << " Us1:( " <<  v1[0] << ","<< v1[1] << ","<< v1[2] << ")" << "    "       
                 << " Us2:( " <<  v2[0] << ","<< v2[1] << ","<< v2[2] << ")" << "    "
                 << " sample: " <<  sampleValue << std::endl;
          }
       
       
     
        // 0 - alphaAlpha
        sampleValue =  alpha1 * alpha2;
       
        insertSample(sampleValue, x1_,1); //pushes the sample into the container
        if(input().verbose())
        {
            std::cout << " at cell: " << celli << std::endl
                 << " alpha1: " <<  alpha1 << "  "
                 << " alpha2: " <<  alpha2 << "  "
                 << " sample: " <<  sampleValue << std::endl;
        }      
        
       }
      }      
                   
    
  }  //end loop over all samples


}

//---------------------------------------------------------------------------------------------

void TwoPointCorr::twoPointsCorrFluct()
{
  int id=dataStorage().fVFid(VFtoSample_[0]);
  if(id==-1) error().throw_error_one("OperationSampling::TwoPointsCorr()",0,"ERROR: Vector field not registered!!"); 
  
  int aId;
  aId=dataStorage().fVFid(averagedField_);
  if(aId==-1) error().throw_error_one("OperationSampling::TwoPointsCorr()",0,"ERROR: averaged Vector field not registered!!"); 
  

  int totcells_;
  double samplePos_[3];
  int celli;
 
  int probeID;
  double v1[3];
  double v2[3];
  double alpha1, alpha2;
  double sampleValue;
  
 
  int cellID;

  double x1_;
  
  if(lagrangian_) dataStorage().setProbes(probesName_);
  
  alpha_=dataStorage().fSFid(alphaName_);
  if(alpha_==-1) error().throw_error_one("OperationSampling::TwoPointsCorr()",0,"ERROR: phase fraction field not registered!!"); 
 
  if(lagrangian_)
   totcells_=dataStorage().numOfParticles();
  else
   totcells_=mesh().NofCells();
   
  

  int j=component_;
    
  for(int rID=0; rID<sampleCount_; rID++)
  {
      x1_= double(rID)*sampleDelta_;
     
      
      for (int cell=0;cell<totcells_;cell++)
      { 
       
        
       if(lagrangian_)
       {
       //calculate CellID
       
       samplePos_[0] = dataStorage().getParticle(cell)->getpos()[0];
       samplePos_[1] = dataStorage().getParticle(cell)->getpos()[1];
       samplePos_[2] = dataStorage().getParticle(cell)->getpos()[2];
       
       cellID = selectorContainer().selectCell(   samplePos_[0],
                                                  samplePos_[1],
                                                  samplePos_[2]
                                              );
       }
       else
       {
        cellID=cell; 
      
        samplePos_[0] = (mesh().CellCentre(cell)[0]);
        samplePos_[1] = (mesh().CellCentre(cell)[1]);
        samplePos_[2] = (mesh().CellCentre(cell)[2]);    
        
        if(dataStorage().SF(alpha_,cellID) > (1-1e-05) ) //no particle in there
         continue;
         
       }
               
       if(cellID>=mesh().NofCells()) error().throw_error_one("OperationSampling::TwoPointsCorr()",0,"ERROR: Problems with cellID!!"); 
      
       //Location for sampling
         
       samplePos_[j] += x1_;

       probeID= selectorContainer().selectCell(   samplePos_[0],
                                                  samplePos_[1],
                                                  samplePos_[2]
                                              );
       
       
       if(probeID>=mesh().NofCells()) 
         error().throw_error_one("OperationSampling::TwoPointsCorr()",0,"ERROR: Problems with probeID!!"); 
       
       if(probeID == -1) continue;      //TODO implement parallel TpointCorr!!!
     
     if(lagrangian_)
       alpha1 = *(dataStorage().getParticle(cell)->filteredScalar(alpha_));
     else
       alpha1 =dataStorage().SF(alpha_,cellID);
       
       alpha2 =dataStorage().SF(alpha_,probeID); 
    
    
        //pull out vector from storage
       
       for (int i=0;i<3;i++)
       {
        if(lagrangian_)
         v1[i]= dataStorage().getParticle(cell)->filteredVector(id)[i] - dataStorage().getParticle(cell)->filteredVector(aId)[i];
        else
         v1[i]= dataStorage().VF(id,i,cellID) - dataStorage().VF(aId,i,cellID);
        
        v2[i]= dataStorage().VF(id,i,probeID) - dataStorage().VF(aId,i,probeID);
       }
       //  cout << "\n" << alpha1 << "\n";
       //  cout << "\n" << alpha2 << "\n";
     
       // 1 - alphaAlphaUsUs 
       sampleValue = alpha1 * alpha2 * (   v1[0] * v2[0] 
                                         + v1[1] * v2[1] 
                                         + v1[2] * v2[2]
                                       );
       if( alpha2 < (1-1e-05) )
       {
          insertSample(sampleValue, x1_); //pushes the sample into the container
          if(input().verbose())
          {
            cout << "\n"
                 << "alpha1: "<< alpha1 << "   alpha2: " << alpha2 << "    " 
                 << " Us1:( " <<  v1[0] << ","<< v1[1] << ","<< v1[2] << ")" << "    "       
                 << " Us2:( " <<  v2[0] << ","<< v2[1] << ","<< v2[2] << ")" << "    "
                 << " sample: " <<  sampleValue << endl;
          }
       
       
     
        // 0 - alphaAlpha
        sampleValue =  alpha1 * alpha2;
       
        insertSample(sampleValue, x1_,1); //pushes the sample into the container
        if(input().verbose())
        {
            cout << " at cell: " << celli << endl
                 << " alpha1: " <<  alpha1 << "  "
                 << " alpha2: " <<  alpha2 << "  "
                 << " sample: " <<  sampleValue << endl;
        }      
        
       }
      }      
                   
    
  }  //end loop over all samples

}

