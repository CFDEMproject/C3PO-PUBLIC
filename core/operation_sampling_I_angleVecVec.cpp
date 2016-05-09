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
#include "operation_sampling.h"
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
#include "operation_sampling_I_angleVecVec.h"

#define PI 3.14159265359

using namespace C3PO_NS;

AngleVecVec::AngleVecVec(c3po *ptr,const char *name) 
: 
OperationSampling(ptr,name),
totalForce_(true)
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /
AngleVecVec::~AngleVecVec()
{
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /

void AngleVecVec::begin_of_step()
{
  (this->*run)();
}

//-----------------------------------------------------------------------------

void AngleVecVec::process_input(QJsonObject jsonObj)
{
 std::string sampleVF(" ");
 std::string sampleSF(" ");
 std::string markers(" ");
 
 if(jsonObj["VecField1"].isNull())
  error().throw_error_one("operation_sampling.cpp",
                                0,
                                "ERROR: You must specify a valid 'VecField1' for angleVecVec sampling \n");
 if(jsonObj["VecField2"].isNull() && jsonObj["particleForce"].isNull() )
  error().throw_error_one("operation_sampling.cpp",
                                0,
                                "ERROR: You must specify a valid 'VecField2' or a valid 'particleForce' for angleVecVec sampling \n");
 if(!jsonObj["VecField2"].isNull() && !jsonObj["particleForce"].isNull() )
  error().throw_error_one("operation_sampling.cpp",
                                0,
                                "ERROR: You can not specify both 'VecField2' and 'particleForce' fields in angleVecVec sampling \n");
 if(!jsonObj["particleForce"].isNull()) lagrangian_=true;

 QString fN=jsonObj["VecField1"].toString();
 sampleVF.assign(fN.toUtf8().constData());
 registerInputFields(sampleVF,sampleSF,markers);

 //for multisampling
 int ssize_=VFtoSample_.size();
 createSampleVectors(ssize_,0,0);
 
 if(!lagrangian_)
 {
  fN=jsonObj["VecField2"].toString();
  sampleVF.assign(fN.toUtf8().constData());
  registerInputFields(sampleVF,sampleSF,markers);
  run=&AngleVecVec::angleVecVec;
 }
 else
 {
  fN=jsonObj["particleForce"].toString();
  if(fN.compare("TotalForce")==0) totalForce_=true;
  else if(fN.compare("ForceModel")==0) totalForce_=false;
  else error().throw_error_one("operation_sampling.cpp",
                                0,
                                "ERROR: You have to enter 'TotalForce' or 'ForceModel' in 'particleForce'\n");
  if(!totalForce_)
  {
  
  if(jsonObj["forceIndex"].isNull() )
   error().throw_error_one("operation_sampling.cpp",
                                0,
                                "ERROR: You have to specify the force index \n");
  forceIndex_=jsonObj["forceIndex"].toInt();
  run=&AngleVecVec::angleVecForce;
                                
  }
   
 }
 
}
//-------------------------------------------------------------------------------
void AngleVecVec::angleVecVec()
{
 int id_2; 
 int markerpos_=VFtoSample_.size()-1;
 id_2=dataStorage().fVFid(VFtoSample_[markerpos_]);
 if(id_2==-1)
  error().throw_error_one("OperationSampling::angleVecVec()",0,"ERROR: Vector field 2 not registered!!");  

 for(int samp_=0; samp_<markerpos_;samp_++)
 {
  int id_1=dataStorage().fVFid(VFtoSample_[samp_]);
  if(id_1==-1)
   error().throw_error_one("OperationSampling::angleVecVec()",0,"ERROR: Vector field 1 not registered!!"); 
 
  double U[3];
  double V[3];

  double angle_;
  double Umod_;

 

//if(save2Bin_)
  //operationContainer().bin(operationID_)->automaticBinCount(3);

  int NofCells_=mesh().NofCells();
 
  for(int cellID=0; cellID<NofCells_; cellID++)
  {
 
 
   for(int j=0;j<3;j++) 
   {
    U[j]=dataStorage().VF(id_1,j,cellID);
    V[j]=dataStorage().VF(id_2,j,cellID); 
   }
  
   Umod_= ( U[0]*U[0] + U[1]*U[1] +U[2]*U[2] );
  
   if(std::fabs(Umod_)<1e-15) return;
  
   angle_ = acos(  std::fabs( ( V[0]*U[0] + V[1]*U[1] + V[2]*U[2])
  
                / sqrt( ( V[0]*V[0] + V[1]*V[1] +V[2]*V[2] ) * Umod_ ) )  
                  
               );
   Umod_=sqrt(Umod_);
   angle_ = angle_ * 180 / PI;  
   if(angle_>180) angle_-= 2*(180 - angle_);                                                                                                                 
 
  
     insertSample(Umod_,angle_,samp_);  
    

  }

 }
 
  flushToBin();
  flushToFile();
  clearSamples();
}
//----------------------------------------------------------------------------------------------------
void AngleVecVec::angleVecForce()
{

 for(unsigned int samp_=0; samp_<VFtoSample_.size();samp_++)
 {
  int id_1=dataStorage().fVFid(VFtoSample_[samp_]);
  if(id_1==-1) error().throw_error_one("OperationSampling::angleVecVec()",0,"ERROR: Vector field 1 not registered!!"); 

  double U[3];
  double V[3];
  double force_[3];
  double angle_;
  double Umod_;

 

//if(save2Bin_)
  //operationContainer().bin(operationID_)->automaticBinCount(3);
  
  if(lagrangian_) dataStorage().setProbes(probesName_);

  int NofParticles_=dataStorage().numOfParticles();
 
  for(int par_=0; par_<NofParticles_; par_++)
  {
       
   
   if(forceIndex_>=dataStorage().getParticle(par_)->getNofForces()) error().throw_error_one("OperationSampling::angleVecVec()",0,"ERROR: No force model corresponding to this forceIndex!!"); 
       
   
 
   for(int j=0;j<3;j++) 
   {
    U[j]=dataStorage().getParticle(par_)->filteredVector(id_1)[j];
    V[j]=dataStorage().getParticle(par_)->getvel()[j];
  
    if(totalForce_) force_[j]=dataStorage().getParticle(par_)->getTotalForce(j);
    else  force_[j]=dataStorage().getParticle(par_)->getforce(forceIndex_)[j];
  
   }
  
   Umod_=( (U[0]-V[0])*(U[0]-V[0]) +(U[1]-V[1])*(U[1]-V[1]) +(U[2]-V[2])*(U[2]-V[2]) );
  
   if(std::fabs(Umod_)<1e-15) return;
  
   angle_ = acos(  std::fabs( ( force_[0]*(U[0]-V[0]) + force_[1]*(U[1]-V[1]) + force_[2]*(U[2]-V[2]))
  
                / sqrt( ( force_[0]*force_[0] + force_[1]*force_[1] +force_[2]*force_[2] ) *  Umod_ )   )   
                  
               );
  
   angle_ = angle_ * 180 / PI;  
   if(angle_>180) angle_-= 2*(180 - angle_);                                                                                                                 

   Umod_=sqrt(Umod_);
   insertSample(Umod_,angle_,samp_);  
  }

 }
 
  flushToBin();
  flushToFile();
  clearSamples();
}

