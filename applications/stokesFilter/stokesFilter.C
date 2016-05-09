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

#include "c3po_OF_interface.H"
#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "treeBoundBox.H"
#include "ListOps.H"
#include "boundBox.H"
#include "direction.H"
#include "pointField.H"
#include "faceList.H"
#include "topoSetSource.H"
#include "treeBoundBoxList.H"
#include "addToRunTimeSelectionTable.H"
#include "dictionary.H"
#include "meshSearch.H"
#include "mpi.h"
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include "error.H" 

#include <ctime>


using namespace Foam;
using namespace C3PO_NS;

#define PI 3.14159265359

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    std::clock_t start;
    start=std::clock();
#   include "setRootCase.H"
#   include "createTime.H"
    
      //MPI initialization
    int flag;
    MPI_Initialized(&flag);
    if (!flag) 
    {
        int argc = 0;
        char **argv = NULL;
        MPI_Init(&argc,&argv);  
    } 
    int proc_name;
    MPI_Comm_rank(MPI_COMM_WORLD,&proc_name);
    
        
   //Read the mesh
    Info<< "Create mesh, no clear-out\n" << endl;
    Foam::fvMesh mesh_
    (
        Foam::IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            Foam::IOobject::MUST_READ
        )
    );
    
    Foam::volVectorField U
    (
        Foam::IOobject
        (
            "U",
            runTime.timeName(),
            mesh_,
            Foam::IOobject::MUST_READ
        ),
        mesh_    
    );
    
    Foam::volScalarField voidfraction
    (
        Foam::IOobject
        (
            "voidfraction",
            runTime.timeName(),
            mesh_,
            Foam::IOobject::MUST_READ
        ),
        mesh_    
    );
    
    double radius_ = 0.5;
    double position_[3];
    double velocity_[3];
    double force_[3];
    
    for(int i=0;i<3;i++)
    {
     position_[i]=-0.0;
     velocity_[i]=0.0;
     force_[i]=0.0;    
    }
   
    
   //  for (int i=0; i<6;i++)
    // std::cout << "\n" << U[i].component(0) << " "<< U[i].component(1) << " "<< U[i].component(2) << " proc: "<< proc_name << " cell: "<< i <<"\n";
   
  
   forAll(mesh_.cells(), celli )
   {
  
   
   //Calculate spherical symmetric coordinates
   
    double x = mesh_.C()[celli].component(0);
    double y = mesh_.C()[celli].component(1);
    double z = mesh_.C()[celli].component(2);
   
    //Info << "\n I am here!";
     
    double r2 = x*x+y*y+z*z;
    double r = std::sqrt(r2);
    
    if (r > radius_) 
    {
     double cosTheta = x/r;
     double tanPhi_;
     double phi;
    if(abs(z)<1e-03)
    {
     if (y<0) phi=-PI/2;
     else phi=PI/2;
    }
    else
    {
     tanPhi_= y/z;
     phi = std::atan( tanPhi_ );     
    }
    
     double theta = std::acos( cosTheta );
     double sinTheta = std::sin(theta);
    //Calculate spherical symmetric components
   
     double U_r =        U[celli].component(0) //the internal field has to be the field at infinite  
                       * cosTheta
                       * (
                             1
                          +  radius_*radius_*radius_ / ( 2*r2*r  )
                          -  3*radius_ / (2*r)
                         );
     
     double U_theta =  -U[celli].component(0) //the internal field has to be the field at infinite  
                       * sinTheta
                       * (
                             1
                          -  radius_*radius_*radius_ / ( 4*r2*r  )
                          -  3*radius_ / (4*r)
                         );
    
    //Calculate Ux, Uy and Uz
   
    U[celli].component(0) = U_r * cosTheta - U_theta * sinTheta;
   
   if(z>0)
   {
     U[celli].component(1) = ( U_r * sinTheta + U_theta * cosTheta )*std::sin(phi);
     U[celli].component(2) =  ( U_r * sinTheta + U_theta * cosTheta )*std::cos(phi);
   }   
   else 
   { U[celli].component(1) = -( U_r * sinTheta + U_theta * cosTheta )*std::sin(phi);
     U[celli].component(2) =  -( U_r * sinTheta + U_theta * cosTheta )*std::cos(phi);
  }
  } 
   }
  // U.correctBoundaryConditions();
   U.write();
   
   
    //Create C3PO
    c3poOFInterface *myC3PO= new c3poOFInterface(mesh_,MPI_COMM_WORLD);
    myC3PO->checkMyMesh();
    
   
   
    std::vector<double*> forceVec_;
    forceVec_.push_back(&force_[0]);
    
    myC3PO->registerParticle(   "particleCenter",
                                0,
                                radius_,
                                &position_[0],
                                &velocity_[0],
                                &forceVec_
                              );
    
    myC3PO->run(); 

    MPI_Barrier(MPI_COMM_WORLD);
    delete myC3PO;
    Pout << "End of application." << endl;
    
    std::clock_t end;
    end=std::clock()-start;
    
    
    
    Pout << "\n the run took: " << (double(end))/CLOCKS_PER_SEC << " s\n";
    
 

}

