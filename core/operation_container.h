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

#ifndef C3PO_OPERATION_CONTAINER_H
#define C3PO_OPERATION_CONTAINER_H

#include "stdio.h"
#include "c3po_base.h"
#include "c3po_base_interface.h"
#include "operation_base.h"
#include <map>
#include <string>
#include "qjson_includes.h"


namespace C3PO_NS
{

class OperationContainer : public c3poBase, public c3poBaseInterface
{
    public:

      OperationContainer(c3po *ptr);
      virtual ~OperationContainer();

      void parse_command(int narg,char const* const* arg);

      void filter_begin_of_step(int s) {opFiltering_[s]->begin_of_step();};
      void filter_middle_of_step(int s) {opFiltering_[s]->middle_of_step();};
      void filter_end_of_step(int s) {opFiltering_[s]->end_of_step();};
      
      void sample_begin_of_step(int s) {opSampling_[s]->begin_of_step();};
      void sample_middle_of_step(int s) {opSampling_[s]->middle_of_step();};
      void sample_end_of_step(int s) {opSampling_[s]->end_of_step();};
      
      void bin_begin_of_step(int s) {opBinning_[s]->begin_of_step();};
      void bin_middle_of_step(int s){opBinning_[s]->middle_of_step();};
      void bin_end_of_step(int s){opBinning_[s]->end_of_step();};
      
      int operationCount() const {return    opFiltering_.size()
                                          + opSampling_.size()
                                          + opBinning_.size();
                                 };

      int filterCount() const {return opFiltering_.size(); };
      int sampleCount() const {return opSampling_.size(); };
      int binCount() const    {return opBinning_.size(); };

      OperationBase* filter(int idx) {return opFiltering_[idx];}
      OperationBase* sample(int idx) {return opSampling_[idx];}
      OperationBase* bin(int idx)    {return opBinning_[idx];}
      void setiOp(int i)             {iOp=i;};
   

      OperationBase* operation(int i)   const ;//return global operation
 
     private:

      typedef OperationBase *(*OperationCreator)(c3po*, char *name);
      std::map<std::string,OperationCreator> *operation_map_;

      template <typename T> static OperationBase *operation_creator(c3po *ptr, char *name);

      std::vector<OperationBase*> opFiltering_;
      std::vector<OperationBase*> opSampling_;
      std::vector<OperationBase*> opBinning_;

      QJsonDocument *jsonDoc_;
      int iOp;
      
};

} //end c3po_NS

#endif
