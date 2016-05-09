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


#include "operation_container.h"
#include "input.h"
#include "style_operation.h"
#include <stdlib.h>
#include <sys/stat.h>
#include "comm.h"
#include "error.h"

#ifdef H5_LIB
#include "h5_c3po.h"
using namespace H5_C3PO_NS;
#endif

using namespace C3PO_NS;

/* ----------------------------------------------------------------------
   OperationContainer Constructor
------------------------------------------------------------------------- */

OperationContainer::OperationContainer(c3po *ptr) 
: 
    c3poBase(ptr),
    operation_map_(new std::map<std::string,OperationCreator>()),
    jsonDoc_(NULL)
  
{
  // fill map with all operations listed in style_operation.h

#define OPERATION_CLASS
#define OperationStyle(key,Class) \
  (*operation_map_)[#key] = &operation_creator<Class>;
#include "style_operation.h"
#undef OperationStyle
#undef OPERATION_CLASS

}

/* ---------------------------------------------------------------------- */

OperationContainer::~OperationContainer()
{
    delete operation_map_;

}

/* ----------------------------------------------------------------------
   one instance per operation in style_operation.h
------------------------------------------------------------------------- */

template <typename T>
OperationBase *OperationContainer::operation_creator(c3po *ptr, char *name)
{
  return new T(ptr,name);
}

/* ----------------------------------------------------------------------
   settings
------------------------------------------------------------------------- */

void OperationContainer::parse_command(int narg,char const* const* arg)
{
    
    int m = strlen(arg[0]) + 1;   
    char *OperationType = new char[m];
    strcpy(OperationType,arg[0]);
    
    int n = strlen(arg[1]) + 1;   
    char *operationName = new char[n];
    strcpy(operationName,arg[1]);
   
      jsonDoc_ = input().mainDoc();
        
        if( jsonDoc_->isNull() )
             error().throw_error_one(FLERR,"QJsonDocument is invalid. Check! \n");
        
        if(jsonDoc_->object()[arg[1]].isNull())
            {
             printf("\nERROR: can not find the operation named %s in c3po.json!\n",arg[1]);
             error().throw_error_one(FLERR,"\nFATAL ERROR\n");
            }
        QJsonObject  operationData = jsonDoc_->object()[arg[1]].toObject();
        
    
        
    if (operation_map_->find(OperationType) != operation_map_->end())
    {
     
        OperationCreator operation_creator = (*operation_map_)[OperationType];
         

        if(strstr(OperationType, "filtering") != NULL)
        {
            opFiltering_.push_back(operation_creator(c3po_ptr(), operationName ) );
            printf("...a filtering operation of type %s is registered with ID %lu\n",
                    OperationType,
                    opFiltering_.size()-1);
            opFiltering_[opFiltering_.size()-1]->init(narg, arg);       //init from input arguments
            opFiltering_[opFiltering_.size()-1]->init(operationData);   //init from JSON file         
        }
        else if(strstr(OperationType, "sampling") != NULL)
        {
            opSampling_.push_back(operation_creator(c3po_ptr(), operationName ) );
            printf("...a sampling operation of type %s is registered with ID %lu\n",
                    OperationType,
                    opSampling_.size()-1);
            opSampling_[opSampling_.size()-1]->init(narg, arg);       //init from input arguments
            opSampling_[opSampling_.size()-1]->init(operationData);   //init from JSON file
        }
        else if(strstr(OperationType, "binning") != NULL)
        {
            opBinning_.push_back(operation_creator(c3po_ptr(), operationName ) );
            printf("...a binning operation is registered with ID %lu\n",
                   opBinning_.size()-1);
            opBinning_[opBinning_.size()-1]->init(narg, arg);       //init from input arguments
            opBinning_[opBinning_.size()-1]->init(operationData);   //init from JSON file
        }
        else
        printf("FAIL: OperationContainer PARSING: operation type not properly set. Use 'filtering', 'sampling', or 'binning' as the type. \n");
    }
    else
        error().throw_error_one(FLERR,"OperationContainer PARSING: operation name not found\n");
        
    

}


// ----------------------------------------------------------------------

OperationBase* OperationContainer::operation(int i) const
{
    int opFiltering_size=opFiltering_.size();
    int opSampling_size=opSampling_.size();
    
    if(i<opFiltering_size)
        return opFiltering_[i];

    else if(i < (opFiltering_size + opSampling_size) )
        return opSampling_[i-opFiltering_size];

    else if(i < operationCount() )
        return opBinning_[i-opFiltering_size-opSampling_size];

    else
        return NULL;

}

