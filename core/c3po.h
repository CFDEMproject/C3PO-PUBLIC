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
	Core class for CPPPO. Provides all the functions needed to communicate with the
	interfaces and holds all tha main classes.
-----------------------------------------------------------------------------------*/


#ifndef C3PO_C3PO_H
#define C3PO_C3PO_H

#include "stdio.h"
#include "mpi.h"
#include <vector>
#include <string>
#include <iostream>
namespace C3PO_NS
{

class c3po
{
 friend class c3poBase;

 public:

  c3po(int narg, char **args, MPI_Comm communicator);
  ~c3po();

  void init();
  void report();
  
 private:

  void help();

  // class to register all members of this class

  class c3poBaseInterfaceVector *c3poBaseInterfaceVector_;

  // ***** CODE FRAMEWORK **************

  // Core Classes
  class Timer               *timer_;
  class Control             *control_;    // controller for stand-alone usage
  class Comm                *comm_;       // inter-processor communication
  class Input               *input_;      // input to simulation - input script processing
  class Output              *output_;     // simulation output
  class Error               *error_;      // error handling
  class OperationContainer  *operationContainer_; // container for filtering/sampling/binning operations
  class SelectorContainer   *selectorContainer_;  // container for selection of cell and particle IDs
  class DataStorage         *dataStorageInternal_;   // storage for particle data, hosts several containers
  class c3poMesh            *mesh_;             // mesh data
  class multiphaseFlowBasic *basicMultiphaseQty_;    // class holding reference quantities and functions (e.g., for drag)

  // ***** END CODE FRAMEWORK **************
    
  public:

    //Misc
    
    //Get number of selectors
    int selectorCount() const;
  
    //Top level run functions
    
    //set cell information
    void   setCells(int * number_of_cells, int * global_number_of_cells, double * cellsize) const;
 
    //run filter id
    void   runFilters(int id) const;

    //run sampling 
    void   runSampling() const;

    //run binning
    void   runBinning() const;
    
    //destroy all operations
    void   flush() const;  //general handle to flush/clear/clean all operation containers

   
   //sampling from interface
      //get number of samples for sampling operation id
       int    sampleCount(int id) ;     
      //get sample delta for sampling operation id
       double sampleDelta(int id) ; 

    //General settings
    bool   verbose() const;

    //mesh getter functions
    double meshCheckTolerance() const;
    double meshFilterWidthTolerance() const;
    bool   meshVerbose() const;
    
    //functions for registering filtered fields 
    void registerVF(std::string,double*,double*,double*,int sp=1);
    void registerSF(std::string,double*);
    
    //functions for registering base fields
    void GlobalVF(std::string , double*,double*,double*, int sp=1);
    void GlobalSF(std::string, double* );
    
    //functions for clearing fields
    void resetFields();
    void resetGlobalFields();
    
    std::string getVFnames(int);
    std::string getSFnames(int);
    
    //get number of fields
    int getVFnamesNumber();
    int getSFnamesNumber();
    
    int getVFnumber();
    int getSFnumber();
    
    //get domain info from json
    double getMaxDomainfromJson(int i) const;
    double getMinDomainfromJson(int i) const;
    
    double cellSizefromJson(int i) const;
    
    
    //add particle from interface
    void addParticle( std::string groupName_,                   //Name of the particle group
                      double   rp,                              //Particle radius    
                      double* pos,                              //Particle position    
                      double* vel,std::vector< double*>* force, //Particle forces
                      std::vector<double>* scalars =NULL,       //Particle scalar fileds (optional)
                      double* torque = NULL                     //Particle torque (optional)
                    );
                    
    //delete particles
    void deleteParticles();
    
    //get name of filtering operation
    const char* getOpFilterName(int id); 
    
    //set simulation time
    void setTime(std::string t);
    
    //get number of operations
    int getOpFiltNum();
    int getOpSampNum();
    int getOpBinNum();
    
    //covert field to global
    void addFieldToConvert(std::string name_);
    
    //refresh fields
    void refreshFields();
    
    //get filter name
    std::string getFilterName(int id);
    
    //get number of filters
    int numberOfFilters();
    
    //get names of fields to filter
    std::vector<std::string> vectorFTF(int id);
    std::vector<std::string> scalarFTF(int id);

    std::vector<std::string> vectorFTFVariance(int id);
    std::vector<std::string> scalarFTFVariance(int id);
    
    //mesh functions
    void registerDomainInfo(double maxDomain[3],double minDomain[3],double maxDomainGlobal[3],double minDomainGlobal[3]);
    void registerCells(double* Vvector, double* posVector, int NumberOfCells, int meshXYZspacing);    
    void setNofCells(const int x) const;
    
    //additional calls
    void preRunOperations() const;
    void postRunOperations() const;
    
    //csv interface functions
    int getDomainDecomposition() const;
    void csv_columnsToRead( int *columns) const;   
    
    //set if mesh is structured
    void checkIJK(bool struct_)   const;
    
    //write to disk
    bool writeFields() const;
    
    //clear mesh data
    void clearMesh() const; 
    
    //true if CFDEM particles are registered
    bool registerCFDEMparticles() const; 
    
    //get names for gradient calculation
    std::vector<std::string> getGradientScalarList() const;
    std::vector<std::string> getGradientVectorList() const;
    std::vector<std::string> getShearRateList() const;
    
    
};

} //end c3po_NS

#endif

