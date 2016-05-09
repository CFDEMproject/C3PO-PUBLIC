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


#include "multiphase_FlowBasic.h"
#include "error.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


using namespace C3PO_NS;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
// Construct from scratch
multiphaseFlowBasic::multiphaseFlowBasic
(
    c3po *ptr
)
:
    c3poBase(ptr),
    verbose(true)
{
    settling.u        = 0.0;
    settling.Re       = 0.0;
    settling.Tref     = 0.0;
    settling.Lref     = 0.0;
    settling.Lchar    = 0.0;
    settling.Lchar2   = 0.0;

    multiphaseSystem.dp         = 0.0;
    multiphaseSystem.rhoP       = 0.0; 
    multiphaseSystem.g          = 0.0;
    multiphaseSystem.etaFluid   = 0.0;
    multiphaseSystem.rhoFluid   = 0.0;    
    multiphaseSystem.dragLaw    = 0.0;  
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //
multiphaseFlowBasic::~multiphaseFlowBasic()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void 
multiphaseFlowBasic::setupSettling( double dP, 
                                    double rhoP, 
                                    double g, 
                                    double etaFluid, 
                                    double rhoFluid, 
                                    int    dragLaw
                                  )
{

    //save multiphase system parameters
    multiphaseSystem.dp         = dP ;
    multiphaseSystem.rhoP       = rhoP; 
    multiphaseSystem.g          = g;
    multiphaseSystem.etaFluid   = etaFluid;
    multiphaseSystem.rhoFluid   = rhoFluid;    
    multiphaseSystem.dragLaw    = dragLaw;  

    //Calculate the settling velocity of an isolated particle
    uSettlingSphereSuspension
                        ( 
                           dP, 
                           rhoP, 
                           g, 
                           etaFluid, 
                           rhoFluid, 
                           0, 
                           dragLaw,
                           settling.u,
                           settling.Re
                        );
    settling.Tref = settling.u / g;
    settling.Lref = settling.u*settling.u / g;
    settling.FrP  = settling.u*settling.u / g / dP;
    settling.Lchar = settling.Lref * std::pow(settling.FrP, -0.5);
    settling.Lchar2= settling.Lref * std::pow(settling.FrP, -0.666666666667);

    if(verbose)
    {           
        printf("\n *** Settling properties have been set up ***\n");
        printf("  dp:     %g [m] \n", multiphaseSystem.dp);
//        printf("  u:      " << settling.u << tab << "[m/s]" << endl;
//        printf("  Re:     " << settling.Re << endl;
//        printf("  Tref:   " << settling.Tref << tab << "[s]" << endl;
//        printf("  Lref:   " << settling.Lref << tab << "[m]" << endl;
//        printf("  FrP:    " << settling.FrP << endl;
//        printf("  Lchar:  " << settling.Lchar << tab << "[m]" <<  endl;
//        printf("  Lchar2: " << settling.Lchar2 << tab << "[m]" <<  endl << endl;
    }
    return;
}; //end function

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void 
multiphaseFlowBasic::uSettlingSphereSuspension
                        ( 
                           double dP, 
                           double rhoP, 
                           double g, 
                           double etaFluid, 
                           double rhoFluid, 
                           double phiP, 
                           int    dragLaw,
                           double &u,
                           double &Re
                          )
{
//uSettlingSphere - Calculates the Settling Velocity of a sphere in a
//homogeneous suspension
//   uSettlingSphereSuspension( dP, rhoP, g, etaFluid, rhoFluid, phiP,
//   dragLaw )
// -------------------------------------
//   INPUT:
//   dP       particle diameter
//   rhoP     particle density
//   g        gravitational acceleration
//   etaFluid dyn. viscosity of fluid
//   rhoFluid density of the fluid
//   phiP     particle volume fraction
//   dragLaw  indicator for drag law
//               0  Wen-Yu
//               1  Beetstra et al.
// -------------------------------------
//   OUTPUT:
//   u        settling velocity
//   Re       Reynolds number


double uInit  = dP * dP * (rhoP-rhoFluid) *g / (18 * etaFluid); //init with stokes settling velocity
double ReInit = fabs(uInit) * dP * rhoFluid / etaFluid;

double relError = 1.0;
double uOld(0.0);
uOld = uInit;
double cDStokes;
double F;

u     = uInit;
Re    = ReInit;

while(relError > 1.0e-6)
{
    
    // Calculate the drag coefficient
    cDStokes = 24/ Re; //Stokes drag coefficient
    if(dragLaw==0)      // Wen-Yu Drag law
        F = F_WenYu(    phiP, Re );
    else if(dragLaw==1)
        F = F_Beetstra( phiP, Re );
    else if(dragLaw==2)
        F = F_KochHill( phiP, Re );
    else
        F = 1.;
    
    u    = sqrt( 4./3. * fabs(rhoP-rhoFluid) * g * dP 
            / ( cDStokes * (1.-phiP) * F * rhoFluid  ) 
             );
             
    Re  = (1.-phiP) * fabs(u) * dP * rhoFluid  / etaFluid;
    
    relError = fabs(uOld - u);
    uOld = u;

}

 return;

}; //end function

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
double
multiphaseFlowBasic::F_Beetstra
          ( 
            double phiP, 
            double Re
          )
{
    double F;
    F    = 10. * phiP / (1.-phiP) / (1.-phiP) 
           + (1.-phiP)*(1.-phiP) * (1+1.5*sqrt(phiP)) 
           +  0.413 * Re / 24. / (1.-phiP) / (1.-phiP)  
            *( 1./(1.-phiP) + 3.0 * phiP * (1-phiP) + 8.4 * std::pow(Re,-0.343) )
            /( 1. + std::pow(10.,3.*phiP) * std::pow(Re,-(1.+4.*phiP)/2.) );

    return F;
}; //end function

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
double
multiphaseFlowBasic::F_WenYu
          ( 
            double phiP, 
            double Re
          )
{

    double F;
	double voidfraction = 1. - phiP;

    F   = 24. * ( 1. + 0.15*std::pow((voidfraction*Re), 0.687)) / (voidfraction*Re);


    return F;
}; //end function

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
double
multiphaseFlowBasic::F_KochHill
          ( 
            double phiP, 
            double Re
          )
{

    double F;
    double F0, F3;
    double voidfraction = 1. - phiP;
        
    //F0
    if(phiP < 0.4)
    {
              F0 = ( 1.+3.*sqrt((phiP)/2)+135./64.*phiP*log(phiP+1.e-16) 
                        +16.14*phiP 
                   ) / 
                   ( 1+0.681*phiP-8.48*phiP*phiP 
                              +8.16*phiP*phiP*phiP 
                   );
    }
    else
    {
              F0 = 10.*phiP/(voidfraction*voidfraction*voidfraction);
    }
        
    //F3
    F3 = 0.0673 + 0.212*phiP + 0.0232 / std::pow(voidfraction,5);
    
    F = voidfraction*(F0 + F3 * Re);

    return F;
}; //end function

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
double
multiphaseFlowBasic::interphaseEC( double particleVolumeFraction, 
                                double gasVelocity, double solidsVelocity
                              )
{

    double slipVelocity = gasVelocity - solidsVelocity;
	double Re = std::abs(slipVelocity) * multiphaseSystem.dp * multiphaseSystem.rhoFluid / multiphaseSystem.etaFluid;
	double voidfraction = 1. - particleVolumeFraction;
    double K = 1;
	

    //use parameters in struct "multiphaseSystem". 
    //TODO This struct must be filled with the function "setupSettling"

    //this is the drag force per m³ suspension
    if(multiphaseSystem.dragLaw == F_WENYU)
    {
     double Cd = F_WenYu(particleVolumeFraction, Re);
     K = (3./4) * Cd * ( particleVolumeFraction * voidfraction * multiphaseSystem.rhoFluid * std::abs(slipVelocity)) * std::pow(voidfraction,-2.65)/ multiphaseSystem.dp;
    }
    else if(multiphaseSystem.dragLaw  == F_BEETSTRA )
    {
     K = 5;
    }
    else if(multiphaseSystem.dragLaw  == F_KOCHHILL )
    {
    }
    else
    {
     error().throw_error_one(FLERR, "\nMultiphaseFlowBasic : drag law does not exist!\nValid drag laws are:\nBeetstra=0\nWenYu=1\nKochHill=2\n");
    }
	 //TODO compute drag force here

    return K;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

 
