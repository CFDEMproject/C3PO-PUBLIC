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
#include "operation_filtering_kernel.h"
#include "selector_container.h"
#include "string.h"
#include "output.h"
#include "input.h"
#include "error.h"
#include "memory_ns.h"
using namespace C3PO_NS;

/* ----------------------------------------------------------------------
   ModelBase Constructor
------------------------------------------------------------------------- */

FilteringKernel::FilteringKernel(c3po *ptr,const char *name) 
: 
OperationFiltering(ptr,name),
weight_fields(0),
inv_weight_fields(0)
{
}

FilteringKernel::~FilteringKernel()
{
 
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void FilteringKernel::process_input(QJsonObject jsonObj)
{
 
 //Check if "kernelSettings" exists
 if(jsonObj["kernelSettings"].isNull())
   error().throw_error_one(FLERR,"You must specify the 'kernelSettings' field for every filtering operation! \n");
   
 //Get kernelSettings object
 QJsonObject kernelSet_ = jsonObj["kernelSettings"].toObject();
 
 
 //Check and register weight fields
 if(!kernelSet_["weightFields"].isNull())
 {
  if(input().verbose())
   output().write_screen_one("\nPerforming weighted average\n");

  
  QString qsS=kernelSet_["weightFields"].toString();
  std::string bigS=qsS.toUtf8().constData();
   
  for (unsigned int it=0; it<bigS.size();it++)
  {
        
        if (bigS[it] != ' ' )
        {  
          
          std::string name;
         
          while (bigS[it] != ' ')
            {
               name.push_back(bigS[it]);
               it++;  
               if (it==bigS.size()) break;           
              
            }  
        //  if(name.empty()==false)  
        // { 
          int index = input().getSFnumber(name);
          if (index==-1) error().throw_error_all("operation_filtering.cpp",0, "Weight fields should be registered as scalar fields in C3PO first");
          
          weightSFName_.push_back(name);
          
         //} 
        }
  }
  
  weight_fields = weightSFName_.size();
  
 }
 
 //Check and register inverted weight fields
 if(!kernelSet_["invertedWeightFields"].isNull())
 {
  if(input().verbose())
   output().write_screen_one("\nPerforming weighted average with inverse\n");

  
  QString qsS=kernelSet_["invertedWeightFields"].toString();
  std::string bigS=qsS.toUtf8().constData();
   
  for (unsigned int it=0; it<bigS.size();it++)
  {
        
        if (bigS[it] != ' ' )
        {  
          
          std::string name;
         
          while (bigS[it] != ' ')
            {
               name.push_back(bigS[it]);
               it++;  
               if (it==bigS.size()) break;           
              
            }  
        //  if(name.empty()==false)  
        // { 
          int index = input().getSFnumber(name);
          if (index==-1) error().throw_error_all("operation_filtering.cpp",0, "Weight fields should be registered as scalar fields in C3PO first");
          { 
           inv_weightSFName_.push_back(name);
          }
         //} 
        }
  }
  
  inv_weight_fields = filterSFName_.size();

 
 }
  
}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void FilteringKernel::check_additional_fields()
{

  weights_id.clear();
  
  for(int i=0;i<weight_fields;i++)
  {
    weights_id.push_back(dataStorage().fSFid(weightSFName_[i]));
    if (weights_id[i]==-1)
        error().throw_error_all("OperationFiltering::begin_of_step()",0,"ERROR: The weight scalar field is NOT registered!");
  }
  
  for(int i=0;i<inv_weight_fields;i++)
  {
    weights_id.push_back(dataStorage().fSFid(inv_weightSFName_[i]));
    if (weights_id[weight_fields+i]==-1)
        error().throw_error_all("OperationFiltering::begin_of_step()",0,"ERROR: The weight scalar field is NOT registered!");
  }
  
  //Create array to hold partial weigths
  //This array is rather expensive in terms of memory
  //But it is necessary for every calculation
  //C3PO_MEMORY_NS::destroy(filterWeights_);
  
  if(!par_)
  {
   filterWeights_= C3PO_MEMORY_NS::create(filterWeights_, mesh().NofCells());
   setVectorToZero(filterWeights_,mesh().NofCells());
  }


}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
inline double FilteringKernel:: calculateWeights( int Id)
{

 //initial weight is the cell volume
 double wi =  *(mesh().CellVol(Id));
 
 //multiply by weights 
 for(int i=0;i<weight_fields;i++)
  wi = wi * dataStorage().SF(weights_id[i],Id);

 //multiply by inverse weights
 for(int i=0;i<inv_weight_fields;i++)
  wi = wi * ( 1.0 - dataStorage().SF(weights_id[weight_fields+i],Id));

 //Return weight
 return wi; 
}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

