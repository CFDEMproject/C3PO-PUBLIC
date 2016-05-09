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
	Core class for multiphase flow calculations and functions. Provides all the 
        reference values and functions needed to normalize key results.
-----------------------------------------------------------------------------------*/
//Contributing Author
//    Stefan Radl (TU Graz, 2015)

#ifndef multiphaseFlowBasic_H
#define multiphaseFlowBasic_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "c3po_base.h"
#include "c3po_base_interface.h"
#include <cmath>

enum{ F_BEETSTRA,
      F_WENYU,
      F_KOCHHILL
     };

namespace C3PO_NS
{

/*---------------------------------------------------------------------------*\
                           Class multiphaseFlowBasic Declaration
\*---------------------------------------------------------------------------*/
class multiphaseFlowBasic : public c3poBase, public c3poBaseInterface
{ 
public:
    bool    verbose;

protected:

    // Protected data

public:

    // Constructor
    multiphaseFlowBasic(c3po *ptr);

    // Destructor
    ~multiphaseFlowBasic();
    
    //variables
    struct settlingParams
    {   
        double  u;              //the settling velocity of an isolated particle
        double  Re;             //the Re of an isolated particle
        double  Tref;           //u/g, the reference time scale
        double  Lref;           //u^2/g, the reference length scale
        double  Lchar;          //u^2/g*FrP^-0.5, a characteristic length scale
        double  Lchar2;          //u^2/g*FrP^-0.5, a characteristic length scale
        double  FrP;             //u^2/g/dp, the particle Froude number
    };

	

    struct systemParams
    {   
        double  dp;             //the settling velocity of an isolated particle
        double  rhoP;           //the Re of an isolated particle
        double  g;              //u/g, the reference time scale
        double  etaFluid;       //u^2/g, the reference length scale
        double  rhoFluid;       //u^2/g*FrP^-0.5, a characteristic length scale
        double  dragLaw;        //u^2/g*FrP^-0.5, a characteristic length scale
    };
    
    settlingParams settling;
    systemParams   multiphaseSystem;
    
    // Member Functions
    double
    F_Beetstra
              ( 
            double, 
            double
              );

    double
    F_WenYu
              ( 
            double, 
            double
              );
              
    double              
    F_KochHill
              ( 
            double, 
            double
              );

    void 
    uSettlingSphereSuspension( double dP, 
                               double rhoP, 
                               double g, 
                               double etaFluid, 
                               double rhoFluid, 
                               double phiP, 
                               int    dragLaw,
                               double &u,
                               double &Re
                               );
                               
    void 
    setupSettling(             double dP, 
                               double rhoP, 
                               double g, 
                               double etaFluid, 
                               double rhoFluid, 
                               int    dragLaw
                 );

    double
    interphaseEC(        double particleVolumeFraction, 
                         double gasVelocity,
			 double solidsVelocity
                );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace C3PO_NS

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

