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
    This class contains those functions needed to make CPPPO to communicate data with
    OpenFOAM. 
-----------------------------------------------------------------------------------*/


#ifndef c3poOFInterface_H
#define c3poOFInterface_H

#include "dictionary.H"
#include "coordinateSystem.H"
#include "coordinateSystems.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedTensor.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "fvMesh.H"
#include "c3po.h"
#include "mesh_check.H"
#include "error.h"
#include <mpi.h>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace C3PO_NS
{

using namespace Foam;

class fvMesh;

/*---------------------------------------------------------------------------*\
Class c3poOFInterface Declaration
\*---------------------------------------------------------------------------*/

class c3poOFInterface
{
    // Private data

private:
        //- Name of this zone, or a regular expression
// Foam::keyType key_;

        //- Reference to the finite volume mesh this zone is part of
        const       Foam::fvMesh& mesh_;
        c3po*       myC3po_;
        MPI_Comm    communicator_;
        

        std::vector< volVectorField*>  vectorFieldsIn_;
        std::vector<const char *>     nameVectorFieldsIn_;

        std::vector< volScalarField*>  scalarFieldsIn_;
        std::vector<const char *>     nameScalarFieldsIn_;
        
        std::vector<volVectorField*>  vectorFieldsC3PO_;
        std::vector<std::string>      nameVectorFieldsC3PO_;
         
        std::vector<volScalarField*>  scalarFieldsC3PO_;
        std::vector<std::string>      nameScalarFieldsC3PO_;
        
        
        std::vector<int>                           particleID_;
        std::vector<double*>                       parVel_;
        std::vector<double>                        pard_;
        std::vector< std::vector<double*> >        parForceTmp_;
        std::vector< std::vector<double> >         parScalarTmp_;
        std::vector<double*>                       parPos_;
        std::vector<double*>                       parTorque_;
        
        
        bool debug_;	
        bool           isPeriodic_[3];
        int me_OF_, nprocs_OF_;
        void processParticles();
        
        meshCheck* mCheck_;  
        
        std::string particlesName_;      
        
        void createGradients();

      

public:

    // Constructors
    //- Construct from components
    c3poOFInterface(const Foam::fvMesh&, MPI_Comm communicator);

    //- Destructor
    ~c3poOFInterface();

    // Member Functions
    void registerVectorField( volVectorField*, const char*);
    void registerScalarField( volScalarField*, const char*);

    int  registeredScalarFields() const {return scalarFieldsIn_.size();};
    int  registeredVectorFields() const {return vectorFieldsIn_.size();};

    //will run all necessary checks on the mesh
    
    void run(); //will run filters, samples, and binning
    void runIB(); //will run filters, samples, and binning for IB solvers (the only difference is in particle force calculation))
    void runC3po();
    void runFilter(int id);
    void runSampling();
    void runBinning();
    
    void checkMyMesh() const {mCheck_->checkMyMesh();};
    void readBC() const {mCheck_->readBC();};
    
    void registerParticle(std::string groupName,
                          int id_,
                          double m,
                          double* pos_,
                          double* vel_,
                          std::vector< double* > * force_,
                          std::vector< double> * scalars_=NULL,
                          double* torque = NULL);
    void registerParticleIB(std::string groupName,
                            int id_,
                            double m,
                            double* pos_,
                            double* vel_,
                            std::vector< double* > force_,
                            std::vector< double> * scalars_=NULL,
                            double* torque = NULL); //parallel handling of immerse boundary forces
    void deleteParticles();

    void registerAllFields();
	
	void registerC3POfields();
	void deleteC3POfields();
	void deleteFields();
	
	
	void writeFields();
	
	void createFilteredFields(int filtId_);
	
	

    bool registerCFDEMparticles() {return myC3po_->registerCFDEMparticles();};
    
    

   // Access

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
