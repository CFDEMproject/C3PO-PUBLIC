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
    
    Foam::volScalarField alpha
    (
        Foam::IOobject
        (
            "alpha",
            runTime.timeName(),
            mesh_,
            Foam::IOobject::MUST_READ
        ),
        mesh_    
    );
    
     Foam::volScalarField p
    (
        Foam::IOobject
        (
            "p",
            runTime.timeName(),
            mesh_,
            Foam::IOobject::MUST_READ
        ),
        mesh_    
    );
    
    //Create harmonics
    
    double PI = 3.14159265359;
    double L = mesh_.bounds().max()[0] - mesh_.bounds().min()[0]; 
    //number of harmonics
    int nmax=1;
    for(int n=1; n<nmax+1; n++)
    {
     forAll(p,cellI)
     {
      
      //wave number
      double k=2*n*PI/(L/2.0);
      
      double x = mesh_.C()[cellI].component(0);
      
      p[cellI] =p[cellI] + (1.0/(2*n))*std::sin(k*x);
      
      U[cellI][0] = U[cellI][0] + (1.0/(2*n))*std::sin(k*x);
      U[cellI][1] = U[cellI][1] + (1.0/(2*n))*std::sin(k*x);
      U[cellI][2] = U[cellI][2] + (1.0/(2*n))*std::sin(k*x);
    
    
     } 
    }    
    p.write();
    U.write();
   alpha.write();
    //  for (int i=0; i<6;i++)
    // std::cout << "\n" << U[i].component(0) << " "<< U[i].component(1) << " "<< U[i].component(2) << " proc: "<< proc_name << " cell: "<< i <<"\n";
    
    
    //Create C3PO
    c3poOFInterface *myC3PO= new c3poOFInterface(mesh_,MPI_COMM_WORLD);
    myC3PO->checkMyMesh();
    myC3PO->run();

//runTime.write();
  //  MPI_Barrier(MPI_COMM_WORLD);
    delete myC3PO;
    Pout << "End of application." << endl;
    
    std::clock_t end;
    end=std::clock()-start;
    
    Pout << "\n the run took: " << (double(end))/CLOCKS_PER_SEC << " s\n";
    
 

}

