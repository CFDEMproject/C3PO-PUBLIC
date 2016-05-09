/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    pisoFoam

Description
    Transient solver for incompressible flow.

    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

//Include CPPPO
#include "c3po_OF_interface.H"

#include "OFversion.H"
#if defined(version30)
    #include "turbulentTransportModel.H"
    #include "pisoControl.H"
#else
    #include "turbulenceModel.H"
#endif


#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "fvIOoptionList.H"
#include "dynamicFvMesh.H" //dyM
#include "cellSet.H"
#include "mpi.h"




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    #if defined(version30)
      pisoControl piso(mesh);
      #include "createTimeControls.H"
    #endif

    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "createFvOptions.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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





 # include "c3po_modifications_1.H"



    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #if defined(version30)
              #include "readTimeControls.H"
              #include "CourantNo.H"
              #include "setDeltaT.H"
          #else
              #include "readPISOControls.H"
              #include "CourantNo.H"
          #endif



        // Pressure-velocity PISO corrector
        {
            // Momentum predictor

            fvVectorMatrix UEqn
            (
                fvm::ddt(U)
              + fvm::div(phi, U)
              + turbulence->divDevReff(U)
              - fvOptions(U)
            );

            UEqn.relax();

            #if defined(version30)
                   if (piso.momentumPredictor())
             #else
                   if (momentumPredictor)
             #endif
            {
                solve(UEqn == -fvc::grad(p));
            }

            // --- PISO loop

            #if defined(version30)
                  while (piso.correct())
              #else
                  int nCorrSoph = nCorr + 5 * pow((1-particleCloud.dataExchangeM().timeStepFraction()),1);
                  for (int corr=0; corr<nCorrSoph; corr++)
              #endif
            {
                volScalarField rAU(1.0/UEqn.A());

                volVectorField HbyA("HbyA", U);
                HbyA = rAU*UEqn.H();
                surfaceScalarField phiHbyA
                (
                    "phiHbyA",
                    (fvc::interpolate(HbyA) & mesh.Sf())
                  + fvc::interpolate(rAU)*fvc::ddtCorr(U, phi)
                );

                adjustPhi(phiHbyA, U, p);

                // Non-orthogonal pressure corrector loop
                #if defined(version30)
                  while (piso.correctNonOrthogonal())
                #else
                  for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                #endif
                {
                    // Pressure corrector

                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rAU, p) == fvc::div(phiHbyA)
                    );

                    pEqn.setReference(pRefCell, pRefValue);

                    #if defined(version30)
                     pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                     #else
                     if( corr == nCorr-1 && nonOrth == nNonOrthCorr )
                       pEqn.solve(mesh.solver("pFinal"));
                     else
                      pEqn.solve();


                     #endif

                }

                #include "continuityErrs.H"

                U = HbyA - rAU*fvc::grad(p);
                U.correctBoundaryConditions();
            }
        }

        turbulence->correct();

//In this example CPPPO is called when OF writes data to disk

     if(runTime.write())
     {
     # include "c3po_modifications_2.H"

        runTime.write();
     }
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

 # include "c3po_modifications_3.H"

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
