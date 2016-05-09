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
#ifdef OPERATION_CLASS

OperationStyle(samplingAngleVecVec,AngleVecVec)

#else

#ifndef C3PO_OPERATION_SAMPLING_I_ANGLEVECVEC_H
#define C3PO_OPERATION_SAMPLING_I_ANGLEVECVEC_H

#include "operation_sampling.h"
#include "output.h"
#include "qjson_includes.h"


namespace C3PO_NS
{
 class AngleVecVec : public OperationSampling
 {
  public:
  
  virtual void process_input(QJsonObject jsonObj);
 
  void begin_of_step();
  void middle_of_step() {};
  void end_of_step() {};
  
  AngleVecVec(c3po *ptr,const char *name);
  ~AngleVecVec();
    
  private:
  
  void angleVecVec();
  void angleVecForce();
  void (AngleVecVec::*run)();
  
  bool totalForce_;
  int forceIndex_;
  
 };


}
#endif
#endif
