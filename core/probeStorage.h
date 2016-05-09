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
/*-----------------------------------------------------------------------------------
Description
     Class for handling groups of particles or probes
-----------------------------------------------------------------------------------*/

#ifndef C3PO_probeStorage_H
#define C3PO_probeStorage_H


#include "c3po_base_accessible.h"
#include "comm.h"
#include <vector>
#include "qjson_includes.h"
#include "particle.h"
#include <string>

namespace C3PO_NS
{
 class probeStorage : public c3poBaseAccessible
 {
  private:
  
  std::string                           groupName_; //identifies the probe group
  
  std::vector<std::string>              filterNames_; //name of the filters that use this group
  
  mutable std::vector<Particle*>        particles_;  //collection of particles
  
  mutable int MaxNofPar_;
        
  mutable std::vector<double>        parPos_; //Holds particle positions for manually registered particles;
 
  mutable int*                       nParticlesProc_;
  
  mutable bool                       readParticlesFromJson_;
  
  public:
 
 
  probeStorage(std::string groupName, std::vector<std::string> filterNames, bool readParticlesFromJson, c3po *ptr);
  ~probeStorage();
 
  void addParticle(double m, double* pos, double* vel, std::vector< double* >* force, std::vector<double>* scalars = NULL, double* torque = NULL);
  
  void gatherParticleData() const;
      
  void readParticles() const;
      
  void removeGhostParticles() const;
  
      
  
  Particle* getParticle(int i) {return particles_[i];};
  
  int numOfParticles() {return particles_.size();};
      
  int MaxNumOfParticles() const { return MaxNofPar_;};
            
  int getNofParticlesProc(int p) const {return  nParticlesProc_[p];};

  
  void writeParticles(std::string fileName_);
  void writeParticleFields(std::string fileName_);
  
  void deleteParticles();
  
  void setParPos(std::vector<double> parPos) const {parPos_=parPos;};
  
  std::string name() {return groupName_;};
  
  bool runProbes(std::string filterName);
  
  void addVector();
  void addScalar();
  
  void cleanProbeFields();
   
 };
}



#endif
