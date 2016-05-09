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
//     Description:
//   sampling class for particles and probes which allows to relate filtered 
//   fields and particle quantities
//
/*-----------------------------------------------------------------------------*/
#ifdef OPERATION_CLASS

OperationStyle(samplingParticles,SamplingParticles)

#else


#ifndef C3PO_OPERATION_SAMPLING_I_PARTICLES_H
#define C3PO_OPERATION_SAMPLING_I_PARTICLES_H

#include "operation_sampling.h"
#include "output.h"
#include "qjson_includes.h"


namespace C3PO_NS
{
 class SamplingParticles : public OperationSampling
 {
  public:
  
  virtual void process_input(QJsonObject jsonObj);
 
  void begin_of_step();

  
  SamplingParticles(c3po *ptr,const char *name);
  ~SamplingParticles();
    
  private:
  
  void particles();
  void checkParticleFields();
  void (SamplingParticles::*run)();
  
  void checkReJson(QJsonObject jsonObj);
  void checkShJson(QJsonObject jsonObj);
  void checkSauterJson(QJsonObject jsonObj);
  
  typedef double(SamplingParticles::*FunctionPointer)(int id);
  
  std::vector< FunctionPointer >   getSampleValue;
  std::vector< FunctionPointer >   getMarkerValue;
  
  double Sh(int id);
  double SauterD(int id);
  double Re(int id);
  double parForce(int id);
  
  void   checkRegisteredFields(std::string name_);
  
  bool    selectRadius_;
  double  maxRad_;
  double  minRad_;                  
  
  int*      VecIds_;
  int*      ScalIds_;
  int*      MarkIds_;
  
  //Limits for samples number (bias statistics)
   double               sample_limits[2];
  
  //  -1   MEMBERS FOR REYNOLDS      ( Re = mod(filtU_Re - Vpar) * dpar * (1 - phi_tot)/ nu_ )
           
           //name of the velocity field
           std::string          filtU_Re;
           
           //kinematic viscosity of the fluid phase
           double               nu_;  
  
           //total fraction of particles
           std::string          phi_tot;   
 

 //  -2   MEMBERS FOR SHERWOOD     ( Sh = (Q*Re*Pr) / dp^2*pi*abs(filtX_Sh-Xref_) )
           
           //name of the scalar field
           std::string          filtX_Sh;
           
           //reference scalar field value
           double               Xref_;  
  
           //id of the particle scalar corresponding to the source
           int                  Qid_;  
           
           //Prandtl number
           double               Pr_;  
           
       
 
 //  -3   MEMBERS FOR SAUTER WEIGHT      ( chi_i = ( phi_par/(phi_tot * dpar) ) / sauterd_ )
          //Actually it is the inverse of the Sauter diameter. 
          //In this way, it is easier to calculate the Sauter mean diameter
          
           //fraction of particles with diameter i
           std::vector<std::string>           phi_i;
           
           //list of diameters
           std::vector<double>                  d_i;
           
           //Id of the current particle diameter
           int                                 d_id;
           //total fraction of particles
           //std::string          phi_tot;   
 
 
 };

}

#endif
#endif
