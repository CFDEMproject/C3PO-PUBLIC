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

#ifndef FILTERED_FIELDS_H
#define FILTERED_FIELDS_H


#include <string>


namespace C3PO_NS
{
 class filteredVectorField 
 {
  private:
  
  std::string fieldName_;
  double * reference_[3]; 
  double filter_size_[3];
  int spacing_;

  int type;
  
  
  
  public:
  
  filteredVectorField(std::string,double*,double*,double*,int spacing=1 );
  filteredVectorField(std::string,int);
  ~filteredVectorField();
  
  void setValue(int i,double* x) {reference_[i]=x;};
  void setFilter(int i, double x) {filter_size_[i]=x;};
  
  double** values() {return reference_;};           //Warning: data in colums (second index) has a certain spacing!
  double*  value(int i,int j) {return &(reference_[i][j*spacing_]);};
  double filter(int i) {return filter_size_[i];};
  std::string name() const {return fieldName_;};
  void setName(std::string x) {fieldName_.assign(x);};
  inline int* space() {return &spacing_;};
  
 };
 
  class filteredScalarField
 {
  private:
  
  std::string fieldName_;
  double * reference_;
  double filter_size_[3];
  int type;  

  
  
  public:
  
  filteredScalarField(std::string,double*);
  filteredScalarField(std::string,int);
  ~filteredScalarField();
  
  void setValue(double* x) {reference_=x;};
  void setFilter(int i, double x) {filter_size_[i]=x;};
  
  double* value() {return reference_;};
  double filter(int i) {return filter_size_[i];};
  std::string name() const {return fieldName_;};
  void setName(std::string x) {fieldName_.assign(x);};
  
 };
}
#endif
