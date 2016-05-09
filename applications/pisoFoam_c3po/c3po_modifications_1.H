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
    
Contributions
     Federico Municchi
\*---------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------
Description
	This file has the sole scope to illustrate how to include the CPPPO libraries
	in an openFOAM solver.
-----------------------------------------------------------------------------------*/
//
//  * * * * * * * * * * * * LINK THE LIBRARIES * * * * * * * * * * * * * * * 

//       Remember to include the following preprocessor directives at the beginning of
//       your solver C file (e.g. pisoFoam_c3po.C) after the header:
//                                         #include "c3po_OF_interface.H"

//
//      And to include the necessary libraries in the Make/options file
                                         
// * * * * * * * * * * * * CODE MODIFICATIONS * * * * * * * * * * * * * * *


// These lines are used to create the C3PO object and run a check of the mesh.In case 
// CPPPO finds incompabilities with the input script an error will be thrown.



 Info << "\n Creating C3PO object...\n";
    
// The c3poOFInterface object requires a fvMesh& and a valid MPI communicator as arguments for 
// its constructor
  
    C3PO_NS::c3poOFInterface *myC3PO= new C3PO_NS::c3poOFInterface(mesh,MPI_COMM_WORLD);
   
    
    
// checkMyMesh() tells CPPPO to read the mesh and collect data. In case of multiple calls, 
// the mesh information will be updated. It has to be called at least one time before doing 
// any operation.       
     
      myC3PO->checkMyMesh();


// Particles can be registered in CPPPO at any point of the solver C file. Normally, if there is no change 
// in particle properties (position, velocity, ecc.) they can be registered just one time before the main loop.
// Note that particles are automatically deleted only when the c3poOFinterface is deleted. Particles can be deleted 
// at any point adding the line:
//                               myC3PO->deleteParticles(); 

 
// Every particle requires a number of parameters. CPPPO accept a std::vector<double*> for the force acting on a particle
// in order to be able to distinguish between different forces (drag, buoyancy, ecc.) 
     
     
 /*    int NofPar= 1;
     double radius=0.25;
     double position[3]={3.9,3.9,3.9};
     double velocity[3]={0,0,0};
     double torque[3]={0,0,0};
     
     double force[3]={0,0,0};
     std::vector<double*> forceVec;
     forceVec.push_back(&force[0]);
     
     for(int par=0;par<NofPar;par++)
     {
     
      myC3PO->registerParticle( par,
                                radius,
                                &position[0],
                                &velocity[0],
                                &forceVec,
                                &torque[0]
                              );
     }                          
     
 */
 
//In addition it is possible to manually define particles using the json files in c3po_control
