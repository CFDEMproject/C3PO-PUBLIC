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

#ifndef C3PO_OPERATION_FILTERING_FILTERINGKERNEL_H
#define C3PO_OPERATION_FILTERING_FILTERINGKERNEL_H

#include "operation_filtering.h"
#include <vector>

namespace C3PO_NS
{
 class FilteringKernel : public OperationFiltering
 {
  public:
  
  FilteringKernel(c3po *ptr,const char *_name);
  ~FilteringKernel();

  protected:
  
  //Names of the weight fields
  mutable std::vector<std::string>    weightSFName_;
  
  //number of weight fields
  mutable int                         weight_fields;
  
  //Names of the inverted weight fields
  mutable std::vector<std::string>    inv_weightSFName_;
    
  //number of inverted weight fields
  mutable int                         inv_weight_fields;
  
  //List of Ids for weight fields
  mutable std::vector<int>                   weights_id;
     
  
  //Check json file for Kernel Settings
  void process_input(QJsonObject jsonObj);
  
  //Virtual function for Kernel implementation
  virtual inline double kernelFunction( double* pos1, double* pos2) {return 1.0;};
  
  //calculate weights
  inline double calculateWeights( int Id);
  
  //Check if weights are registerd
  void check_additional_fields();
   
 };
}
#endif
