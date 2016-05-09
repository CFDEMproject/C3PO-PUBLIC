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
///
//  * * * * * * * * * * * * LINK THE LIBRARIES * * * * * * * * * * * * * * * 

//       Remember to include the following preprocessor directives at the beginning of
//       your solver C file (e.g. pisoFoam_c3po.C) after the header:
//                                         #include "c3po_OF_interface.H"

//
//      And to include the necessary libraries in the Make/options file
                                         
// * * * * * * * * * * * * CODE MODIFICATIONS * * * * * * * * * * * * * * *



    Info << "Running C3PO post-processing";
    
// This line will tell CPPPO to process your data accordingly to the input files located in
// the c3po_control directory. Most of the applications will require the following line.
   
     myC3PO->run(); 
        
// In case of a parallel DNS solver, the same particle could be subject to different forces from 
// different processors. In this situation (that occurs in cfdemIBSolver) it is necessary to use a
// different CPPPO function and substituting myC3PO->run() with:
//                                                     myC3PO->runIB(); 

     

