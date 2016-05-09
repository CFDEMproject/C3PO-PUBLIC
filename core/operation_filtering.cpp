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

#include "operation_filtering.h"
#include "operation_container.h"
#include "selector_container.h"
#include "filter_base.h"
#include "string.h"
#include "comm.h"
#include "output.h"
#include "error.h"
#include <fstream>
#include <cmath>
#include <string>
#include <sstream>

#include "timer.h"

#include "memory_ns.h"

#define VERBOSE_OPFILTER false

using namespace C3PO_NS;
using namespace std;

/* ----------------------------------------------------------------------
   ModelBase Constructor
------------------------------------------------------------------------- */

OperationFiltering::OperationFiltering(c3po *ptr,const char *name) 
: 
OperationBase(ptr,name),
probesName_("NULL"),
computeVariance_(false),
varianceHasCrossTerm_(false),
par_(false),
VF_to_filter(0),
SF_to_filter(0),
VF_for_varianceCalc_(0),
SF_for_varianceCalc_(0)
{

    operationID_ = operationContainer().filterCount(); //set ID, this is the 'old' sampleCount!
}

OperationFiltering::~OperationFiltering()
{
 
 
 delete fieldsToAdd_;
 
 C3PO_MEMORY_NS::destroy(fieldsPerProc_);
 C3PO_MEMORY_NS::destroy(tmpData_);
 C3PO_MEMORY_NS::destroy(VTmp_);
 C3PO_MEMORY_NS::destroy(Var_);
 C3PO_MEMORY_NS::destroy(filterWeights_);
}

// * * * * * * * * * * * * * * * * * * * * * * * *
void OperationFiltering::init(QJsonObject jsonObj) 
{

 QJsonObject     mainSettings_;
 QJsonObject varianceSettings_;
 
  
 if(jsonObj["mainSettings"].isNull())
   error().throw_error_one(FLERR,"You must specify the \"mainSettings\" in your filtering operation!!. \n");
 
 mainSettings_=  jsonObj["mainSettings"].toObject();
 
  if(mainSettings_["lagrangian"].isNull())
   error().throw_error_one(FLERR,"You must specify the \"lagrangian\" field!!. \n");
         
  par_=mainSettings_["lagrangian"].toBool();
  
  if(par_)
  {
    if(mainSettings_["probesName"].isNull())
     error().throw_error_one(FLERR,"If you set a filtering operation in \"lagrangian\" mode you should also specify the \"probesName\" field! \n");
     
     QString qsV=mainSettings_["probesName"].toString();
     probesName_ =qsV.toUtf8().constData();
   
  }
  
   if(!jsonObj["varianceSettings"].isNull())
   {
    
    varianceSettings_=jsonObj["varianceSettings"].toObject();
    
    if(!varianceSettings_["computeVariance"].isNull())
     computeVariance_ =  varianceSettings_["computeVariance"].toBool();
   }  
 
   registerFieldsFiltered(mainSettings_); 
  
  if(computeVariance_)
   registerFieldsVariance(varianceSettings_); 
 
 
   totFields_=SF_to_filter + 3*(VF_to_filter) +1;
  
   
   fieldsPerProc_= C3PO_MEMORY_NS::create(fieldsPerProc_, comm().nprocs(), totFields_);
   tmpData_=C3PO_MEMORY_NS::create(tmpData_, comm().nprocs(), totFields_);
 
   if(computeVariance_)            
   {
    totVar_=SF_for_varianceCalc_+(VF_for_varianceCalc_*3);
    VTmp_= C3PO_MEMORY_NS::create(VTmp_, comm().nprocs(), totVar_);
    Var_= C3PO_MEMORY_NS::create(Var_, comm().nprocs(), totVar_);
   
   } 
  
  if(!par_)
  {
   run=&OperationFiltering::divergentAlgorithm;
   end=&OperationFiltering::endDivergentAlgorithm; 
  }
  else
  {
   run=&OperationFiltering::convergentAlgorithm;
   end=&OperationFiltering::endConvergentAlgorithm; 
  }
  
   
  process_input(jsonObj);
  
  
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ 
void OperationFiltering::registerFieldsFiltered(QJsonObject jsonObj) const
{
   QString qsV=jsonObj["VectorfieldsToFilter"].toString();
   std::string big=qsV.toUtf8().constData();
  
  for (unsigned int it=0; it<big.size();it++)
  {
        
        if (big[it] != ' ' )
        {  
          
          std::string name;
         
          while (big[it] != ' ')
            {
               name.push_back(big[it]);
               it++;  
               if (it==big.size()) break;           
              
            }  
       // if(name.empty()==false)  
        // {
          int index = input().getVFnumber(name);
          if (index==-1) error().throw_error_all("operation_filtering.cpp",0, "Vector fields to filter should be registered in C3PO first");
          filterVFName_.push_back(name);
        // }  
        }
  }
  
  VF_to_filter=filterVFName_.size();
  
  QString qsS=jsonObj["ScalarfieldsToFilter"].toString();
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
          if (index==-1) error().throw_error_all("operation_filtering.cpp",0, "Scalar fields to filter should be registered in C3PO first");
          { 
           filterSFName_.push_back(name);
          }
         //} 
        }
  }
  
  SF_to_filter = filterSFName_.size();

}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ 
void OperationFiltering::registerFieldsVariance(QJsonObject jsonObj) const
{
  //Vector Fields
  QString qsV=jsonObj["VectorfieldsForVarianceName1"].toString();
  std::string bigV=qsV.toUtf8().constData();
  for (unsigned int it=0; it<bigV.size();it++)
  {
        if (bigV[it] != ' ' )
        {  
          std::string name;
          while (bigV[it] != ' ')
          {
               name.push_back(bigV[it]);
               it++;  
               if (it==bigV.size()) break;           
          }  
          int index = input().getVFnumber(name);
          if (index==-1) error().throw_error_all("operation_filtering.cpp",0, "Vector fields to filter should be registered in C3PO first");
          filterVFVarianceName_.push_back(name);

        }
  }
  VF_for_varianceCalc_=filterVFVarianceName_.size();

  //Set names for second value in (OPTIONAL) for cross-variance calculation
  qsV=jsonObj["VectorfieldsForVarianceName2"].toString();
  bigV=qsV.toUtf8().constData();
  for (unsigned int it=0; it<bigV.size();it++)
  {
        if (bigV[it] != ' ' )
        {  
          std::string name;
          while (bigV[it] != ' ')
          {
               name.push_back(bigV[it]);
               it++;  
               if (it==bigV.size()) break;           
          }  
          int index = input().getVFnumber(name);
          if (index==-1) error().throw_error_all(FLERR, "Vector fields to filter should be registered in C3PO first");
          filterVFVarianceName2_.push_back(name);
        }
  }
  //check for correct input
  if(filterVFVarianceName2_.size()>0 && (filterVFVarianceName_.size()!=filterVFVarianceName2_.size()) )
    error().throw_error_all(FLERR, "You DID NOT specify enough/too much values for VectorfieldsForVarianceName2");
  
  if(filterVFVarianceName2_.size()==0)
  {
    for(uint iVar = 0; iVar<filterVFVarianceName_.size(); iVar++)
        filterVFVarianceName2_.push_back(filterVFVarianceName_[iVar]);
  }
  else   varianceHasCrossTerm_ = true;

  //Decide to compute off-diagonal elements (and not diagonal elements)
  QJsonArray tmpArray=jsonObj["VectorfieldsForVarianceComputeOffDiagonal"].toArray();
  //check for correct input
  if(tmpArray.size()>0 && ((uint)tmpArray.size()!=filterVFVarianceName_.size()) )
    error().throw_error_all(FLERR, "You DID NOT specify enough/too much bools for VectorfieldsForVarianceComputeOffDiagonal");
  
  if(tmpArray.size()==0)
  {
    for(uint iVar = 0; iVar<filterVFVarianceName_.size(); iVar++)
       filterVFVarianceComputeOffDiagonal_.push_back(false);
  }
  else
  {
    for (int it=0; it<tmpArray.size();it++)
        filterVFVarianceComputeOffDiagonal_.push_back( (tmpArray.at(it)).toBool() );
  }
 
  //Scalar Fields
  QString qsS=jsonObj["ScalarfieldsForVarianceName1"].toString();
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
          int index = input().getSFnumber(name);
          if (index==-1) error().throw_error_all(FLERR, "Scalar fields to filter should be registered in C3PO first");
          filterSFVarianceName_.push_back(name);
          computeVariance_ = true;
        }
  }
  SF_for_varianceCalc_ = filterSFVarianceName_.size();

  //Set names for second value in (OPTIONAL) for cross-variance calculation
  qsS=jsonObj["ScalarfieldsForVarianceName2"].toString();
  bigS=qsS.toUtf8().constData();
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
          int index = input().getSFnumber(name);
          if (index==-1) error().throw_error_all(FLERR, "Scalar fields to filter should be registered in C3PO first");
          filterSFVarianceName2_.push_back(name);
        }
  }
  //check for correct input
  if(filterSFVarianceName2_.size()>0 && (filterSFVarianceName_.size()!=filterSFVarianceName2_.size()) )
    error().throw_error_all(FLERR, "You DID NOT specify enough/too much values for ScalarfieldsForVarianceName2");

  if(filterSFVarianceName2_.size()==0)
    for(uint iVar = 0; iVar<filterSFVarianceName_.size(); iVar++)
        filterSFVarianceName2_.push_back(filterSFVarianceName_[iVar]);
  else
    varianceHasCrossTerm_ = true;



  //MIXED: section for reading vector-scalar mixed variance calculations
  qsS=jsonObj["ScalarfieldsForVectorScalarMixedVariance"].toString();
  bigS=qsS.toUtf8().constData();
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
          int index =  searchScalarFiltName(name);
          if (strcmp(name.c_str(),"off")==0)
            std::cout << "Found 'off' setting, hence ScalarfieldsForVectorScalarMixedVariance will NOT be performed for component " 
                      << filterSFVarianceNameVectorScalarMixed_.size() << ".\n";
          else if(index==-1)
          {
            std::cout << "ERROR: Could not find a 'ScalarfieldsToFilter' with name '" << name << "'! Be sure to add this scalar field to 'ScalarfieldsToFilter' in case you want to evaluate vector-scalar mixed variances.";
            error().throw_error_all(FLERR, "Must abort.");
          }
          filterSFVarianceNameVectorScalarMixed_.push_back(name);
        }
  }
  //check for correct input
  if(filterSFVarianceNameVectorScalarMixed_.size()>0 && (filterVFVarianceName_.size()!=filterSFVarianceNameVectorScalarMixed_.size()) )
    error().throw_error_all(FLERR, "You DID NOT specify enough/too much values for ScalarfieldsForVectorScalarMixedVariance.");

  if(varianceHasCrossTerm_ && filterSFVarianceNameVectorScalarMixed_.size()>0)
    error().throw_error_all(FLERR, "You cannot specify 'VectorfieldsForVarianceName2' and 'ScalarfieldsForVectorScalarMixedVariance'. Please use separate filtering operations in case you want to compute vector-vector cross terms and vector-scalar mixed terms.");

  for(unsigned int iVar = 0; iVar<filterSFVarianceNameVectorScalarMixed_.size(); iVar++)
  {
    int tmpScalInt = searchScalarFiltName(filterSFVarianceNameVectorScalarMixed_[iVar]);
    filterSFVarianceValueVectorScalarMixed_.push_back(tmpScalInt);
    if(tmpScalInt>=0)
        evaluateVarianceVectorScalarMixed_.push_back(true);
    else
        evaluateVarianceVectorScalarMixed_.push_back(false);

    if(filterVFVarianceComputeOffDiagonal_[iVar])
        std::cout << "WARNING: you have set 'VectorfieldsForVarianceComputeOffDiagonal' to true. This is irrelevant for vector-scalar mixed variance calculation.\n";
#if VERBOSE_OPFILTER
    std::cout << "filterSFVarianceNameVectorScalarMixed_[" << iVar << "]: " 
              << filterSFVarianceNameVectorScalarMixed_[iVar] 
              << ": identified the following scalar id " << tmpScalInt
              << ", hence bool for vector-scalar mixed is: " << evaluateVarianceVectorScalarMixed_[iVar] << "\n";
#endif
    
  }
}



/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void OperationFiltering::begin_of_step()
{

 //Filter fields
 for(int i=0;i<VF_to_filter;i++)
  {

   RMAVF_.push_back(dataStorage().fVFid(filterVFName_[i]));
   if (RMAVF_[i]==-1)
    error().throw_error_all("OperationFiltering::begin_of_step()",0,"ERROR: The vector field you want to filter is NOT registered!");
   filterVF_.push_back(dataStorage().VFid(filterVFName_[i]));
    if (filterVF_[i]==-1)
    error().throw_error_all("OperationFiltering::begin_of_step()",0,"ERROR: The vector field you want to filter is NOT registered!!");
   std::string name(dataStorage().fVF(filterVF_[i])->name());
    name.append("_");
    name.append(name_);
    dataStorage().fVF(filterVF_[i])->setName(name);
  }
  
  for(int i=0;i<SF_to_filter;i++)
  { 

    RMASF_.push_back(dataStorage().fSFid(filterSFName_[i]));
    if (RMASF_[i]==-1)
        error().throw_error_all("OperationFiltering::begin_of_step()",0,"ERROR: The scalar field you want to filter is NOT registered!");
    filterSF_.push_back(dataStorage().SFid(filterSFName_[i]));
    if (filterSF_[i]==-1)
        error().throw_error_all("OperationFiltering::begin_of_step()",0,"ERROR: The scalar field you want to filter is NOT registered!!");
    std::string name(dataStorage().fSF(filterSF_[i])->name());
    name.append("_");
    name.append(name_);
    dataStorage().fSF(filterSF_[i])->setName(name);

#if VERBOSE_OPFILTER
    std::cout << "filterSFName_[" << i << "]: " 
              << filterSFName_[i] 
              << ": added the following field to dataStorage: " << dataStorage().fSF(filterSF_[i])->name()
              << " in slot: " << filterSF_[i] << "\n";
#endif

  }

  //Variance Fields - VECTOR Quantities
  for(int j=0;j<VF_for_varianceCalc_;j++)
  {
    int i = VF_to_filter + j;
    filterVFVarianceValueID_.push_back(       searchVectorFiltName(filterVFVarianceName_[j]) );
    filterVFVarianceValueSecondID_.push_back( searchVectorFiltName(filterVFVarianceName2_[j]) );
        
    filterVF_.push_back(dataStorage().VFid(filterVFVarianceName_[j]));
    if (filterVF_[i]==-1)
        error().throw_error_all("OperationFiltering::begin_of_step()",
                                -1,
                                "ERROR: The vector field you want to use for variance calculation is NOT registered in VFid!");
        
    std::string name(dataStorage().fVF(filterVFVarianceName_[j])->name());
    name.append("_");
    name.append(name_);

    //Check corresponding filtered field id to verify registration of this field
    int filtFieldID = dataStorage().VFid(name);
    if (filtFieldID==-1)
    {
        error().throw_error_all("OperationFiltering::begin_of_step()",
                                -1,
                                "ERROR: Could not find filtered VECTOR field for variance calculation. Check you have the above filtered field registerd!");
    }
    name.append("_variance");
    ostringstream Str; Str << j;  name.append(Str.str());
    dataStorage().fVF(filterVF_[i])->setName(name);

#if VERBOSE_OPFILTER
    std::cout << "filterVFVarianceName_[" << j << "]: " 
              << filterVFVarianceName_[j] 
              << ": in dataStorage, value ID: " << filterVFVarianceValueID_[j] 
              << ", valueSecond ID: " << filterVFVarianceValueSecondID_[j]  << "\n"; 
    std::cout << "filterVFVarianceName_[" << j << "]: " 
              << filterVFVarianceName_[j] 
              << ": added the following field to dataStorage: " << dataStorage().fVF(filterVF_[i])->name()
              << " in slot: " << filterVF_[i] << "\n";
#endif
  }


  //Variance Fields - SCALAR Quantities
  for(int j=0;j<SF_for_varianceCalc_;j++)
  {
    int i = SF_to_filter + j;
    filterSFVarianceValueID_.push_back(       searchScalarFiltName(filterSFVarianceName_[j]) );
    filterSFVarianceValueSecondID_.push_back( searchScalarFiltName(filterSFVarianceName2_[j]) );
        
    filterSF_.push_back(dataStorage().SFid(filterSFVarianceName_[j]));
    if (filterSF_[i]==-1)
        error().throw_error_all("OperationFiltering::begin_of_step()",
                                -1,
                                "ERROR: The scalar field you want to use for variance calculation is NOT registered in SFid!");
        
    std::string name(dataStorage().fSF(filterSFVarianceName_[j])->name());
    name.append("_");
    name.append(name_);

    //Check corresponding filtered field id to verify registration of this field
    int filtFieldID = dataStorage().SFid(name);
    if (filtFieldID==-1)
    {
        std::cout << "ERROR: when searching for filtered field " << name << "\n" ;
        error().throw_error_all("OperationFiltering::begin_of_step()",
                                -1,
                                "ERROR: Could not find filtered SCALAR field for variance calculation. Check you have the above filtered field registerd!");
    }
    name.append("_variance");
    ostringstream Str; Str << j;  name.append(Str.str());
    dataStorage().fSF(filterSF_[i])->setName(name);

#if VERBOSE_OPFILTER
    std::cout << "filterSFVarianceName_[" << j << "]: " 
              << filterSFVarianceName_[j] 
              << ": in dataStorage, value ID: " << filterSFVarianceValueID_[j] 
              << ", valueSecond ID: " << filterSFVarianceValueSecondID_[j]  << "\n"; 
    std::cout << "filterSFVarianceName_[" << j << "]: " 
              << filterSFVarianceName_[j] 
              << ": added the following field to dataStorage: " << dataStorage().fSF(filterSF_[i])->name()
              << " in slot: " << filterSF_[i] << "\n";
#endif
  }
  
  check_additional_fields();
  
  
 
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void OperationFiltering::middle_of_step()
{
 (this->*run)();
}
/* * * * ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */ 
void OperationFiltering::end_of_step()
{
 (this->*end)();
}
/* * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * * */ 
void OperationFiltering::setVectorToZero(double* vector, int DIM)
{
 for(int i=0;i<DIM;i++)
  vector[i]=0.0;
}
/* * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * * */ 
void OperationFiltering::setVectorToOne(double* vector, int DIM)
{
 for(int i=0;i<DIM;i++)
  vector[i]=1.0;
}
/* * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * * */ 
bool OperationFiltering::particleBased()
{
 if(par_)
 {
 
   if(probesName_.compare(dataStorage().currentProbes())==0 || probesName_.compare("all")==0) 
    return true;
 
 }
 return false;
}
