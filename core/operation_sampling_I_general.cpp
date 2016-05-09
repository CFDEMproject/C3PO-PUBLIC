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
#include "operation_sampling_I_general.h"


using namespace C3PO_NS;


SamplingGeneral::SamplingGeneral(c3po *ptr,const char *name) 
: 
OperationSampling(ptr,name)
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /
SamplingGeneral::~SamplingGeneral()
{
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /

void SamplingGeneral::begin_of_step()
{
  (this->*run)();
}

//-----------------------------------------------------------------------------

void SamplingGeneral::process_input(QJsonObject jsonObj)
{
 if(lagrangian_) run=&SamplingGeneral::sampleGeneralLagrangian;
 else            run=&SamplingGeneral::sampleGeneralEulerian;
}

/*------------------------------------------------------------------------- */
void SamplingGeneral::sampleGeneralEulerian()
{ 
   int totcells_;
  double sampleValue, marker1;
  int cellID,id;
      //Checks
 
  if(sampleCount_<0)
        totcells_=mesh().NofCells();
  else
        totcells_=sampleCount_;
  std::vector<int> idMarker_;
  for( int mark_=0;mark_<NofMarkers_;mark_++)
  {
    int idMarker1=dataStorage().fSFid(markers_[mark_]);
    if(idMarker1==-1) 
     error().throw_error_one(FLERR,"Marker field not registered!!"); 
   idMarker_.push_back(idMarker1);
  }       
  

 //process Vector Fields
  for(unsigned int samp_=0;samp_<VFtoSample_.size();samp_++)
  {
   id=dataStorage().fVFid(VFtoSample_[samp_]);
   if(id==-1) 
    error().throw_error_one(FLERR,"Vector field not registered!!"); 
 
 
    for (int cell=0;cell<totcells_;cell++)
    { 
      double markerValue[NofMarkers_];
       cellID=cell;    
              
       if(cellID>=mesh().NofCells()) 
            error().throw_error_one(FLERR,"ERROR: cellID out of range!!"); 
       if(selective_)
       {
        bool skip_=false;
        for(int i=0;i<3;i++)
         if((mesh().CellCentre(cell)[i])>maximum_[i] || (mesh().CellCentre(cell)[i])<minimum_[i] )
          skip_=true;
        
        if(skip_) continue;
       }
         
        sampleValue  = dataStorage().VF(id, component_,cellID);
      
        for( int mark_=0;mark_<NofMarkers_;mark_++)  
           markerValue[mark_]=dataStorage().SF(idMarker_[mark_],cellID);

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
   id=dataStorage().fSFid(SFtoSample_[samp_]);
   if(id==-1) 
    error().throw_error_one(FLERR,"ERROR: Scalar field not registered!!"); 
   
    for (int cell=0;cell<totcells_;cell++)
    { 
       double markerValue[NofMarkers_];
       cellID=cell;    
       
       if(selective_)
       {
        bool skip_=false;
        for(int i=0;i<3;i++)
         if((mesh().CellCentre(cell)[i])>maximum_[i] || (mesh().CellCentre(cell)[i])<minimum_[i] )
          skip_=true;
        
        if(skip_) continue;
       }
              
       if(cellID>=mesh().NofCells()) 
            error().throw_error_one(FLERR,"ERROR: cellID out of range!!"); 
         
        sampleValue  = dataStorage().SF(id,cellID);
        
        for( int mark_=0;mark_<NofMarkers_;mark_++)  
            markerValue[mark_]=dataStorage().SF(idMarker_[mark_],cellID);

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
/*------------------------------------------------------------------------- */
void SamplingGeneral::sampleGeneralLagrangian()
{ 
  int totPar_;
  double sampleValue, marker1;
  int id;
   
  //set correct probe
  dataStorage().setProbes(probesName_);
 
 
  if(sampleCount_<0)
        totPar_=dataStorage().numOfParticles();
  else
       totPar_=sampleCount_;
  
  if(totPar_>dataStorage().numOfParticles()) totPar_=dataStorage().numOfParticles();     
  
  std::vector<int> idMarker_;
  for( int mark_=0;mark_<NofMarkers_;mark_++)
  {
    int idMarker1=dataStorage().SFid(markers_[mark_]);
    if(idMarker1==-1) 
     error().throw_error_one(FLERR,"Marker field not registered!!"); 
   idMarker_.push_back(idMarker1);
  }       
  

 //process Vector Fields
  for(unsigned int samp_=0;samp_<VFtoSample_.size();samp_++)
  {
   id=dataStorage().VFid(VFtoSample_[samp_]);
   if(id==-1) 
    error().throw_error_one(FLERR,"Vector field not registered!!"); 
 
 
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
         
        sampleValue  = dataStorage().getParticle(par)->filteredVector(id)[component_];
      
        for( int mark_=0;mark_<NofMarkers_;mark_++)  
           markerValue[mark_]=*dataStorage().getParticle(par)->filteredScalar(idMarker_[mark_]);

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
   id=dataStorage().SFid(SFtoSample_[samp_]);
   if(id==-1) 
    error().throw_error_one(FLERR,"ERROR: Scalar field not registered!!"); 
   
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
              
   
         
        sampleValue  = *dataStorage().getParticle(par)->filteredScalar(id);
        
        for( int mark_=0;mark_<NofMarkers_;mark_++)  
            markerValue[mark_]=*dataStorage().getParticle(par)->filteredScalar(idMarker_[mark_]);

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
