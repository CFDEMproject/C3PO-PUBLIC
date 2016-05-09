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
    Class to parse and store mesh data from CSV.
-----------------------------------------------------------------------------------*/
#ifndef CSV_LAGRANGIAN_C3PO
#define CSV_LAGRANGIAN_C3PO
#include "c3po.h"
#include  <string>
#include  <fstream>
#include  <vector>
#include "mesh_check.h"

namespace C3PO_NS
{
 class CSVlagrangian
 {
  public:
  
  CSVlagrangian(c3po* , CSVmesh*);
  ~CSVlagrangian();
  
  void checkLagrangian() const;
  void deleteParticles() const;
 
  private:
  
  c3po *      C3po_;
  CSVmesh *   mesh_;
  
  mutable int NofPar_;
  
  mutable double tolerance_;
  
  mutable std::vector< double* > particlePos_;
  mutable std::vector< double* > particleVel_;
  
  mutable std::vector< double* > particleTorque_;
  
  mutable std::vector< std::vector< double* > > particleForce_;
  
  mutable std::vector< double > particleRadius_;
     
  void readLagrangian() const;
  void registerC3POparticles() const; 
  
  
 
 };


}

#endif
