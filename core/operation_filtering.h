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

#ifndef C3PO_OPERATION_FILTERING_H
#define C3PO_OPERATION_FILTERING_H

#include "operation_base.h"
#include <vector>

namespace C3PO_NS
{


class OperationFiltering : public OperationBase 
{
  public:

  OperationFiltering(c3po *ptr,const char *_name);
  ~OperationFiltering();

  void begin_of_step();
  void middle_of_step();
  void end_of_step();
  
  void init(QJsonObject jsonObj);
  
  void registerFieldsFiltered(QJsonObject jsonObj) const;
  
  void registerFieldsVariance(QJsonObject jsonObj) const;
      
  virtual void process_input(QJsonObject jsonObj) {};
    
  bool particleBased(); 
  mutable std::string  probesName_;
      
  //The filtered fields
  std::vector<std::string> returnVectorFTFNames() {return filterVFName_;};
  std::vector<std::string> returnScalarFTFNames() {return filterSFName_; };       

  //The variance fields
  std::vector<std::string> returnVectorFTFVarianceNames() {return filterVFVarianceName_;};
  std::vector<std::string> returnScalarFTFVarianceNames() {return filterSFVarianceName_; };       
  
  protected:
  
  void setVectorToZero(double* vector, int DIM);
  void setVectorToOne(double* vector, int DIM);
  
  virtual void check_additional_fields() {};
    
  mutable std::vector<int> filterVF_;
  mutable std::vector<int> filterSF_;
   
  mutable std::vector<int> RMAVF_;
  mutable std::vector<int> RMASF_;
    
  mutable std::vector<std::string> filterVFName_;
  mutable std::vector<std::string> filterSFName_;

  mutable std::vector<std::string> filterVFVarianceName_;
  mutable std::vector<std::string> filterSFVarianceName_;
  mutable std::vector<std::string> filterVFVarianceName2_;
  mutable std::vector<std::string> filterSFVarianceName2_;
  mutable std::vector<int>  filterVFVarianceValueID_;
  mutable std::vector<int>  filterSFVarianceValueID_;
  mutable std::vector<int>  filterVFVarianceValueSecondID_;
  mutable std::vector<int>  filterSFVarianceValueSecondID_;
  
  mutable std::vector<bool> filterVFVarianceComputeOffDiagonal_;
  mutable bool    computeVariance_;
  mutable bool    varianceHasCrossTerm_;
  mutable int totFields_; 
  mutable double** fieldsPerProc_;
  mutable double** tmpData_;
  
  mutable double** VTmp_;
  mutable double** Var_;
  mutable int totVar_;
  
  mutable double* fieldsToAdd_;
  
  std::vector<double*>  valueContainerScalars_;
  std::vector<double*>  valueContainerVectors_;
  
  void (OperationFiltering::*run)();
  void (OperationFiltering::*end)();

  //for vector-scalar variance calculations
  mutable std::vector<bool>        evaluateVarianceVectorScalarMixed_;
  mutable std::vector<std::string> filterSFVarianceNameVectorScalarMixed_;
  mutable std::vector<int>         filterSFVarianceValueVectorScalarMixed_;


  mutable int VF_to_filter;
  mutable int SF_to_filter; 
  mutable bool par_;

  mutable int VF_for_varianceCalc_;
  mutable int SF_for_varianceCalc_;
 
  double *    filterWeights_; 
  

  int searchScalarFiltName(std::string name) const
  {
        for (unsigned int i=0;i<filterSFName_.size();i++)
        if (filterSFName_[i]==name) return i;
        return -1;
  }

  int searchVectorFiltName(std::string name) const
  {
        for (unsigned int i=0;i<filterVFName_.size();i++)
        if (filterVFName_[i]==name) return i;
        return -1;
  }
  
  
  void divergentAlgorithm();
  void endDivergentAlgorithm();
  
  void convergentAlgorithm();
  void endConvergentAlgorithm() {};
  
  virtual inline double kernelFunction( double* pos1, double* pos2) {return 1.0;};
  virtual inline double calculateWeights(int Cell)            {return 1.0;};
  
        
};

} //end c3po_NS

#endif
