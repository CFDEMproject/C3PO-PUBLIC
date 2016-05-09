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
     Main class to manage references to lagrangian and eulerian data. It contains
     std::vectors of filteredField and Particle objects. 
-----------------------------------------------------------------------------------*/


#ifndef C3PO_DATA_STORAGE_H
#define C3PO_DATA_STORAGE_H

#include "stdio.h"
#include "c3po_base.h"
#include "c3po_base_interface.h"
#include <vector>
#include "filtered_fields.h"
#include "qjson_includes.h"
#include "particle.h"
#include "probeStorage.h"
#include "region.h"


using std::vector;

namespace C3PO_NS
{

class DataStorage : public c3poBase, public c3poBaseInterface
{
    public:

      DataStorage(c3po *ptr);
      ~DataStorage();

      
      void init() const;
      
      void addfVF(std::string,double*,double*,double*,int sp=1);
      
      void addfSF(std::string,double*);
      
      void deleteFields();
      void deleteRMA();
      
      inline filteredVectorField* fVF(int i) {return fVF_[i];};
      inline filteredScalarField* fSF(int i) {return fSF_[i];};
      
      filteredVectorField* RMAvF(int i) const {return RMAvF_[i];};
      filteredScalarField* RMAsF(int i) const {return RMAsF_[i];};
      
      filteredVectorField* fVF(std::string);
      filteredScalarField* fSF(std::string);
      
      int fVFid(std::string name)
      {
        for (unsigned int i=0;i<RMAvF_.size();i++)
        if (RMAvF_[i]->name()==name) return i;
      
        return -1;
      }; 
      
       int fSFid(std::string name)
      {
        for (unsigned int i=0;i<RMAsF_.size();i++)
        if (RMAsF_[i]->name()==name) return i;
        
        return -1;
      
      }; 
      
      int VFid(std::string name)
      {
        for (unsigned int i=0;i<fVF_.size();i++)
        if (fVF_[i]->name()==name) return i;
      
        return -1;
      }; 
      
       int SFid(std::string name)
      {
        for (unsigned int i=0;i<fSF_.size();i++)
        if (fSF_[i]->name()==name) return i;
        
        return -1;
      
      };                         
      
      int numberOfVF() {return RMAvF_.size();};
      int numberOfSF() {return RMAsF_.size();};
      
      void addRMAvF(std::string,double*,double*,double*, int sp=1);
      void addRMAsF(std::string,double*);
      
  
      
      inline double  VF(int name, int component, int cell ) { return *(RMAvF_[name]->value(component,cell)); };
      inline double  SF(int name, int cell ) { return RMAsF_[name]->value()[cell]; };
      
      
      int fVFsize()  {return fVF_.size();};
      int fSFsize()   {return fSF_.size();};

      void writeFields(std::string OpName);
    
      
      std::string NameChanges(std::string OpName);
      void setFileName(std::string OpName);
      
     
      
      void setTimeName(std::string t) const {timeName_.assign(t);};
      std::string getTimeName()       const {return timeName_;};
      
      void convertFilterToRMA(std::string name_);
      void addFieldToConvert(std::string field_);
      void refreshRMAfields();
      
      //Probes/Particles functions
      
      
      void addParticle(std::string groupName, double m, double* pos, double* vel, std::vector< double* >* force, std::vector<double>* scalars = NULL, double* torque = NULL);
     
      Particle* getParticle(int i) {return probes_[currentProbe_]->getParticle(i);};
      
      int numOfParticles() {return probes_[currentProbe_]->numOfParticles();};
      
      int MaxNumOfParticles() const { return  probes_[currentProbe_]->MaxNumOfParticles();};
      
      int getNofParticlesProc(int p) const {return probes_[currentProbe_]->getNofParticlesProc( p);};
      
      int NofProbes() const {return probes_.size();};
      
      void setCurrentProbes(int i)  const {currentProbe_=i;};
      
      void setProbes(std::string probeName_) const;
      
      std::string currentProbes()  { return probes_[currentProbe_]->name();};
      
      void writeParticles();
      void writeParticleFields(std::string filtName);
      
      void readParticles() const;
      
      bool runProbes(std::string filterName) {return probes_[currentProbe_]->runProbes(filterName);};
      
      void deleteParticles();
      
      bool useProbes()  {return useProbes_;};
      
      void clearProbeFields() {probes_[currentProbe_]->cleanProbeFields();};

     
     
      void writeRegions(); //Function to open output file, and loop through regions such that each region "donate" data to the output

      void createRegion(std::string regionName); 
      
      void addCellToRegion(int & IDInternal, bool isOnSurface);
      
      void selectRegion(int id) {regionPtr_=regions_[id];};
      
      void deleteRegions();
      
      int NofRegions() {return regions_.size();};
      

    private:

        mutable std::vector<std::string>          fieldsToConvert_;
        
        std::vector<filteredVectorField*>         fVF_;
        std::vector<filteredScalarField*>         fSF_;
        
        std::vector<filteredVectorField*>         RMAvF_;
        std::vector<filteredScalarField*>         RMAsF_;
       
        mutable std::vector<probeStorage*>        probes_;
        
        //List of Regions     
        std::vector<region*>                      regions_;
        
        region*                                 regionPtr_;
    
        
        int                                 buf1;
        double                              buf2;
        int                                 nproc;
        std::string                         filename_;
        std::string                         dirNameField_;
        std::string                         dirNameParticle_;
        std::string                         dirNameRegions_;
        mutable std::string                 timeName_;
        
        mutable int                               currentProbe_;

        int nbody_, nbody_all_;

        mutable bool isAllocated_;
        bool haveFieldDir_;
        bool haveParticleDir_;
        bool haveRegionsDir_;
        
        mutable bool useProbes_;

       
};

} //end c3po_NS

#endif
