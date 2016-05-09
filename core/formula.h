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
\*-----------------------------------------------------------------------------
    Description
    tool for interpretation and execution of user-defined formulas
\*-----------------------------------------------------------------------------*/
#ifndef C3PO_FORMULA_H
#define C3PO_FORMULA_H

#include "stdio.h"
#include <string>
#include <vector>
#include "mpi.h"
#include <iostream>
#include "multiphase_FlowBasic.h"

namespace C3PO_NS
{
 class Formula
 {
  public:
  
  Formula(const char* formula);
  Formula(const char* formula,const char* normalization);

  ~Formula();
 
  void            interpretFormula(int, int) const; //interprets the std::strings in the formula
  std::string     getFormula() const; //just for display
  std::string     getNumerator() const; //just for display
  std::string     getDenominator() const; //just for display

  void            setLimiterNumerator(double min, double max) const    //set value
                  { limitNumerator_[0] = min; limitNumerator_[1] = max; limitNumeratorActive_ = true;};
  double*         getLimitNumerator() const  {return limitNumerator_;}


  void            setLimiterDenominator(double min, double max) const  //set value
                  { limitDenominator_[0] = min; limitDenominator_[1] = max; limitDenominatorActive_ = true;};
  double*         getLimitDenominator() const  {return limitDenominator_;}

  void            setLimiterSymmetry(bool numerator, bool denominator)
                  {limitNumeratorSymmetric_ = numerator; limitDenominatorSymmetric_=denominator;};

  void            evaluate( std::vector< std::vector<double>  > * inputFields_, 
                            std::vector< std::vector<double>* >*  _markers,
                            std::vector< double > * outputField_,
                            std::vector< std::vector<double> >*    newMarkers_
                          );

  void            normalize( std::vector< std::vector<double>  > * inputFields_, 
                             std::vector< std::vector<double>* >*  _markers,
                             std::vector< double > * outputField_ 
                           ); //normalizes the sample

  void            setMultBasicObject(multiphaseFlowBasic* _myBasic ) {basicMultiphaseQty_ = _myBasic;}; //just to hand over pointer to basic multiphase calculator
  
  
  void            setMaxNumberOfFields(int maxMark_, int maxScal_,int maxVec_) const;
  
  void            setAbsoluteValue(bool absoluteNumerator, bool absoluteDenominator) const {
                                                                                            absoluteNumerator_   =   absoluteNumerator;
                                                                                            absoluteDenominator_ = absoluteDenominator;
                                                                                           };
  

  private:
  
  mutable bool              limitNumeratorActive_;
  mutable bool              limitDenominatorActive_;
  mutable double            limitNumerator_[2];
  mutable bool              limitNumeratorSymmetric_;
  mutable double            limitDenominator_[2];
  mutable bool              limitDenominatorSymmetric_;
  mutable bool              normalize_;
  mutable bool              absoluteNumerator_;
  mutable bool              absoluteDenominator_;
  
  mutable int               NofVecSample_ ;
  mutable int               NofScalSample_;
  mutable int               NofMarkSample_;

  mutable std::string       raw_formula_;
  mutable std::string       normalization_;
  
  mutable std::vector<char> formula_;
  mutable std::vector<char> numerator_; 
  mutable std::vector<char> denominator_;   //by default empty
  
  mutable std::vector<int>  tmpIndex_;
  
  mutable multiphaseFlowBasic* basicMultiphaseQty_;
  
  mutable std::vector<int>  normalizationId_;
  mutable std::vector<char> normalizationClass_;
  mutable char              normalizationType_;
  
  

  void ParseFormula() const;
  void CheckFormula() const;
  void CheckNormalization(int, int) const;
  
  void throw_error(std::string function,std::string msg_) const;
  bool testLimits(double& value, double* limits, bool& isSymmetric) const;
  
  double getValueForNormalization(int id,
                                  int index,
                                  std::vector< std::vector<double> > *  _sample,
                                  std::vector< std::vector<double>* >*  _markers 
                                  ) const;
  
 };
}
#endif
