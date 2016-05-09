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
#include "selector_container.h"
#include <fstream>
#include <sys/stat.h>
#include <iomanip>
#include <cmath>
#include "input.h"
#include "output.h"
#include "operation_sampling_I_particles.h"

#ifdef H5_LIB
#include "h5_c3po.h"
#define PI 3.14159265359
using namespace H5_C3PO_NS;
#endif

using namespace C3PO_NS;


SamplingParticles::SamplingParticles(c3po *ptr,const char *name) 
: 
OperationSampling(ptr,name),
selectRadius_(false),
maxRad_(1e32),
minRad_(-1)
{
 sample_limits[0] = -1e32;
 sample_limits[1] =  1e32;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /
SamplingParticles::~SamplingParticles()
{
 delete VecIds_;
 delete ScalIds_;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /

void SamplingParticles::begin_of_step()
{
  dataStorage().setProbes(probesName_);
  checkParticleFields();
  (this->*run)();
}

//-----------------------------------------------------------------------------

void SamplingParticles::process_input(QJsonObject jsonObj)
{
 //Set the pointer to function 
 run=&SamplingParticles::particles;
 VecIds_= new int[VFtoSample_.size()];
 ScalIds_= new int[SFtoSample_.size()];
 MarkIds_= new int[markers_.size()];
 if(!lagrangian_)
    error().throw_error_one(FLERR,"\n \"lagrangian\" should be set to \"true\" for particle sampling operations! \n");

  if(!jsonObj["samplesLimiter"].isNull())
   {
     QJsonArray limitsTemp=jsonObj["samplesLimiter"].toArray();
     sample_limits[0] = limitsTemp.at(0).toDouble();
     sample_limits[1] = limitsTemp.at(1).toDouble();
   }


 if(!jsonObj["maxParticleRadius"].isNull())
 {
  selectRadius_=true;
  maxRad_ =jsonObj["maxParticleRadius"].toDouble();
 
 }
 
 if(!jsonObj["minParticleRadius"].isNull())
 {
  selectRadius_=true;
  minRad_ =jsonObj["minParticleRadius"].toDouble();
 
 }
 
 //Check requirements for special samples/Markers
 for(unsigned int mark_=0;mark_<markers_.size();mark_++)
 {
  
   if(markers_[mark_].compare("Reynolds")==0)
   {
    checkReJson(jsonObj);
    getMarkerValue.push_back(&SamplingParticles::Re);
   }
   else if(markers_[mark_].compare("SauterW")==0)
   {
    checkSauterJson( jsonObj);
    getMarkerValue.push_back(&SamplingParticles::SauterD);
  
   }
   else if(markers_[mark_].compare("Sherwood")==0)
   {
    checkShJson(jsonObj);
    getMarkerValue.push_back(&SamplingParticles::Sh);
   
   }
   else if(markers_[mark_].compare("particleForce")==0)
   {
    getMarkerValue.push_back(&SamplingParticles::parForce);   
   }
   else
    getMarkerValue.push_back(&SamplingParticles::Re); //push back pointer to keep the ordering
  
 }
 
  for(unsigned int samp_=0;samp_<SFtoSample_.size();samp_++)
  {  
   if(SFtoSample_[samp_].compare("Reynolds")==0)
   {
    checkReJson(jsonObj);
    getSampleValue.push_back(&SamplingParticles::Re);
   }
   else if(SFtoSample_[samp_].compare("SauterW")==0)
   {
    checkSauterJson(jsonObj);
    getSampleValue.push_back(&SamplingParticles::SauterD);
  
   }
   else if(SFtoSample_[samp_].compare("Sherwood")==0)
   {
    checkShJson(jsonObj);
    getSampleValue.push_back(&SamplingParticles::Sh);
   
   }
   else if(SFtoSample_[samp_].compare("particleForce")==0)
   {
    getSampleValue.push_back(&SamplingParticles::parForce);   
   }
   else
    getSampleValue.push_back(&SamplingParticles::Re); //push back pointer to keep the ordering
 
  }
}

//-----------------------------------------------------------------------------

void SamplingParticles::particles()
{

 int totPar_;
  double sampleValue, marker1;
  int id;
   

  if(sampleCount_<0)
        totPar_=dataStorage().numOfParticles();
  else
       totPar_=sampleCount_;
  
  if(totPar_>dataStorage().numOfParticles()) totPar_=dataStorage().numOfParticles();          
  

 //process Vector Fields
  for(unsigned int samp_=0;samp_<VFtoSample_.size();samp_++)
  {
  
    id=VecIds_[samp_];
    
    for (int par=0;par<totPar_;par++)
    { 
      double markerValue[NofMarkers_];  
              
       if(selective_)
       {
        bool skip_=false;
        for(int i=0;i<3;i++)
         if(dataStorage().getParticle(par)->getpos()[i] > maximum_[i]
                        ||
           dataStorage().getParticle(par)->getpos()[i]<minimum_[i]
           )
          skip_=true;
        
        if(skip_) continue;
       }
       
       if(selectRadius_)
       {
        if(*dataStorage().getParticle(par)->getradius()<minRad_ || *dataStorage().getParticle(par)->getradius() > maxRad_)
         continue;
       }
         
        sampleValue  = dataStorage().getParticle(par)->filteredVector(id)[component_];
        
        if(sampleValue < sample_limits[0] || sampleValue > sample_limits[1] )
         continue;
       
        for( int mark_=0;mark_<NofMarkers_;mark_++)  
        {
         if(MarkIds_[mark_]<0) markerValue[mark_]=(this->*(getMarkerValue[mark_]))(par);
         else                  markerValue[mark_]=*dataStorage().getParticle(par)->filteredScalar(MarkIds_[mark_]);
        }
        
      
         
        insertSample(sampleValue,&markerValue[0],samp_); //pushes the sample into the container
        if(input().verbose())
        {
            std::cout<< "VectorField:sampleValue: "<< sampleValue;
            for(int mark_=0;mark_<NofMarkers_;mark_++)
            {
             std::cout<< "   marker" << mark_ << ": " << markerValue[mark_];
            }
            std::cout<< std::endl;
        }
    }
   }
   

 //process Scalar Fields
  for(unsigned int samp_=0;samp_<SFtoSample_.size();samp_++)
  {
   id=ScalIds_[samp_];
   
  
    for (int par=0;par<totPar_;par++)
    { 
       double markerValue[NofMarkers_];
       
       if(selective_)
       {
        bool skip_=false;
        for(int i=0;i<3;i++)
        if(dataStorage().getParticle(par)->getpos()[i] > maximum_[i]  ||
           dataStorage().getParticle(par)->getpos()[i]<minimum_[i]
           )   skip_=true;
        
        if(skip_) continue;
       }
       
       if(selectRadius_)
       {
        if(*dataStorage().getParticle(par)->getradius()<minRad_ || *dataStorage().getParticle(par)->getradius() > maxRad_)
         continue;
       }
              
        if(ScalIds_[samp_]<0) sampleValue  = (this->*(getSampleValue[samp_]))(par);
        else               sampleValue  = *dataStorage().getParticle(par)->filteredScalar(id);
        
        if(sampleValue < sample_limits[0] || sampleValue > sample_limits[1] )
         continue;     
        
        for( int mark_=0;mark_<NofMarkers_;mark_++)  
        {
         if(MarkIds_[mark_]<0) markerValue[mark_]=(this->*(getMarkerValue[mark_]))(par);
         else                  markerValue[mark_]=*dataStorage().getParticle(par)->filteredScalar(MarkIds_[mark_]);
        }
        
        insertSample(sampleValue,&markerValue[0],VFtoSample_.size()+samp_); //pushes the sample into the container
        if(input().verbose())
        {
            std::cout<< "ScalarField:sampleValue: "<< sampleValue;
            for(int mark_=0;mark_<NofMarkers_;mark_++)
            {
             std::cout<< "   marker" << mark_ << ": " << markerValue[mark_];
            }
            std::cout<< std::endl;
        }
    }
  }  

} 
//-----------------------------------------------------------------------------

void SamplingParticles::checkParticleFields()
{
 //Analyse input to distinguish between particle properties (e.g radius)
 //and filtered fields
  
  for( int mark_=0;mark_<NofMarkers_;mark_++)
  {
    int idMarker1=dataStorage().SFid(markers_[mark_]);
    if(idMarker1==-1) //If it is not a field, it may be a particle scalar property 
     checkRegisteredFields(markers_[mark_]);
   
   MarkIds_[mark_]=idMarker1;
  }   
  
  for(unsigned int s_=0;s_<SFtoSample_.size();s_++)
  {
    int id=dataStorage().SFid(SFtoSample_[s_]);
    if(id==-1) //If it is not a field, it may be a particle scalar property 
     checkRegisteredFields(SFtoSample_[s_]);
  
   ScalIds_[s_]=id;
  }  
  
 
  for(unsigned int v_=0;v_<VFtoSample_.size();v_++)
  {
    int id=dataStorage().VFid(VFtoSample_[v_]);
    if(id==-1) //If it is not a field, it may be a particle scalar property 
    {

      error().throw_error_one(FLERR,"Marker field not registered!!"); 
    }
   
   
   VecIds_[v_]=id;
  }            
      
  
}
//-----------------------------------------------------------------------------

double SamplingParticles::Re(int par)
{
  int Uid = dataStorage().VFid(filtU_Re);
  
  int phitot_id = dataStorage().SFid(phi_tot);
  
  double phitot = *dataStorage().getParticle(par)->filteredScalar(phitot_id);
  
  double UrelSq =0.0; 
  
  for(int j=0;j<3;j++)
  {
    double rel =  dataStorage().getParticle(par)->filteredVector(Uid)[j] - dataStorage().getParticle(par)->getvel()[j];
    rel *= rel;
    UrelSq += rel;
  }
 
  return sqrt(UrelSq)*(*dataStorage().getParticle(par)->getradius())*(1.0-phitot)*2/nu_;

}
//-----------------------------------------------------------------------------

double SamplingParticles::Sh(int par)
{
  double Re_ = Re(par);
  int Xid = dataStorage().SFid(filtX_Sh);
  
  double deltaT = Xref_ - *dataStorage().getParticle(par)->filteredScalar(Xid)  ;
  
  if( deltaT < 1e-15 && deltaT > -1e-15) return 0.0; //return no exchange
   
  double dimlessArea = (*dataStorage().getParticle(par)->getradius())
                       *
                       (*dataStorage().getParticle(par)->getradius())
                       *PI*4; 
  
  return (dataStorage().getParticle(par)->getscalar(Qid_)*Pr_*Re_) /(dimlessArea*deltaT);

}
//-----------------------------------------------------------------------------

double SamplingParticles::SauterD(int par)
{
  
  int phitot_id = dataStorage().SFid(phi_tot);
  double phitot = *dataStorage().getParticle(par)->filteredScalar(phitot_id);
  double sauter=0.0;
  if( phitot < 1e-08 ) return 0.0; //no particles here!
  
  for(unsigned int i=0; i<phi_i.size();i++)
  {
   int phii_id = dataStorage().SFid(phi_i[i]);
  
   double phii = *dataStorage().getParticle(par)->filteredScalar(phii_id);
   
   sauter += phii/(d_i[i]*phitot);
   
  }
  
   int phii_id = dataStorage().SFid(phi_i[d_id]);
  
   double phii = *dataStorage().getParticle(par)->filteredScalar(phii_id);
  
  if( sauter < 1e-08 ) return 0.0; //weird!

  return ( phii/(d_i[d_id]*phitot) ) / sauter;

}
//-----------------------------------------------------------------------------

double SamplingParticles::parForce(int par)
{
  
  double forceMod = 0.0;
  
  for(int coord=0; coord<3; coord++)
   forceMod += dataStorage().getParticle(par)->getTotalForce(coord) 
               * dataStorage().getParticle(par)->getTotalForce(coord);
  
  forceMod = sqrt(forceMod);
  return forceMod;
}

//-----------------------------------------------------------------------------
void SamplingParticles::checkReJson(QJsonObject jsonObj)
{
 //Check quantities for Reynolds
 
  if(!jsonObj["filteredVelocityForRe"].isNull())
    filtU_Re = jsonObj["filteredVelocityForRe"].toString().toUtf8().constData();
   else
    error().throw_error_one(FLERR,"Velocity field needed for particle Reynolds number!!"); 
  
  
   if(!jsonObj["nu"].isNull())
    nu_ = jsonObj["nu"].toDouble();
   else
    error().throw_error_one(FLERR,"Fluid viscosity needed for particle Reynolds number!!"); 
  
   if(!jsonObj["particleFractionTot"].isNull())
    phi_tot = jsonObj["particleFractionTot"].toString().toUtf8().constData();
   else
    error().throw_error_one(FLERR,"The superficial Reynolds number requires to specify the total particle fraction!!"); 
  
  
   if(nu_<1e-10)
         error().throw_error_one(FLERR,"Fluid kinematic viscosity \"nu\" cannot be negative or < 1e-10!"); 

 
}

//-----------------------------------------------------------------------------
void SamplingParticles::checkShJson(QJsonObject jsonObj)
{
  //Particle Reynolds number is required
    if(filtU_Re.compare("empty")==0)
     error().throw_error_one(FLERR,"You need a Reynolds number marker to calculate the particle Sherwood number!! "); 
   
   checkReJson(jsonObj);
   
    if(!jsonObj["filteredXForSh"].isNull())
    filtX_Sh = jsonObj["filteredXForSh"].toString().toUtf8().constData();
   else
    error().throw_error_one(FLERR,"Scalar field needed for particle Sherwood number!!"); 
  
  
   if(!jsonObj["Xref"].isNull())
    Xref_ = jsonObj["Xref"].toDouble();
   else
    error().throw_error_one(FLERR,"You need to specify a reference X for the particle Sherwood number!!"); 
    
   if(!jsonObj["Pr"].isNull())
    Pr_ = jsonObj["Pr"].toDouble();
   else
    error().throw_error_one(FLERR,"You need to specify the Prandtl number for the particle Sherwood number!!"); 
  
   if(!jsonObj["Qid"].isNull())
    Qid_ = jsonObj["Qid"].toInt();
   else
    error().throw_error_one(FLERR,"You need to specify the id of the particle scalar property corresponding to the interphase heat to calculate the particle Sherwood number!!"); 
  

}

//-----------------------------------------------------------------------------
void SamplingParticles::checkSauterJson(QJsonObject jsonObj)
{
   if(!jsonObj["particleFractionTot"].isNull())
    phi_tot = jsonObj["particleFractionTot"].toString().toUtf8().constData();
   else
    error().throw_error_one(FLERR,"Total particle fraction field necessary for Sauter diameter!!"); 
  
   if(!jsonObj["particleFraction_i"].isNull())
   {
     QJsonArray  fractionList = jsonObj["particleFraction_i"].toArray();
    for(int i=0;i<fractionList.count();i++)
     phi_i.push_back(fractionList.at(i).toString().toUtf8().constData());
    
   }
   else
    error().throw_error_one(FLERR,"Particle fraction field array for diameter D necessary for Sauter diameter!!"); 
    
   if(!jsonObj["particleDiameter_i"].isNull())
   {
     QJsonArray  fractionList = jsonObj["particleDiameter_i"].toArray();
    for(int i=0;i<fractionList.count();i++)
     d_i.push_back(fractionList.at(i).toDouble());
    
   }
   else
    error().throw_error_one(FLERR,"Particle diameter array necessary for Sauter diameter!!"); 
  
   if(!jsonObj["fraction_id"].isNull()) d_id=jsonObj["fraction_id"].toInt();
   else
    error().throw_error_one(FLERR,"Particle identificator necessary for Sauter diameter!!"); 
  
   if(phi_i.size()!= d_i.size() || d_i.size() < d_id+1)
     error().throw_error_one(FLERR,"Inconsistent input for sauter diameter"); 

}
//-----------------------------------------------------------------------------
void SamplingParticles::checkRegisteredFields(std::string name_)
{
  
    if(name_.compare("Reynolds")==0)
     {
      if(dataStorage().VFid(filtU_Re)==-1)
        error().throw_error_one(FLERR,"Velocity field for particle Reynolds not registered!!"); 
      
       if(dataStorage().SFid(phi_tot)==-1)
        error().throw_error_one(FLERR,"Total particle fraction not registered!!"); 
      
      
     }
     else if(name_.compare("SauterW")==0)
     {
      if(dataStorage().SFid(phi_tot)==-1)
        error().throw_error_one(FLERR,"Total particle fraction not registered!!"); 
      
      for(unsigned int i=0;i<phi_i.size();i++)
      if(dataStorage().SFid(phi_i[i])==-1)
        error().throw_error_one(FLERR,"Particle fraction for diameter not registered!!"); 
   
      
      
     }
     else if(name_.compare("Sherwood")==0)
     {
      if(dataStorage().SFid(filtX_Sh)==-1)
        error().throw_error_one(FLERR,"Scalar field for particle Sherwood not registered!!"); 
      
      if(dataStorage().numOfParticles()>0)
       if(Qid_ > dataStorage().getParticle(0)->checkScalars()-1)
         error().throw_error_one(FLERR,"Not enough scalar fields registered in particle!!"); 
      
      
    }
    else if(name_.compare("particleForce")==0)
    {
    }
    else
     error().throw_error_one(FLERR,"Field not registered!!"); 

}

