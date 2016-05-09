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

#include "filtered_fields.h"
#include <iostream>

using namespace C3PO_NS;

filteredVectorField::filteredVectorField(std::string name, double* x, double* y, double* z, int spacing /*=1*/)
:
 fieldName_(name), 
 spacing_(spacing),
 type(0)
{ 

 reference_[0]=x;
 reference_[1]=y;
 reference_[2]=z;

 filter_size_[0]=-1;
 filter_size_[1]=-1;
 filter_size_[2]=-1;

}

filteredVectorField::filteredVectorField(std::string name,int number)
:
 fieldName_(name), 
 spacing_(1),
 type(1)
{ 

 reference_[0]=new double[number];
 reference_[1]=new double[number];
 reference_[2]=new double[number];

 filter_size_[0]=-1;
 filter_size_[1]=-1;
 filter_size_[2]=-1;

}

 
filteredVectorField::~filteredVectorField()
{
  if (type==1)
  {
   for (int i=0;i<3;i++)
    delete reference_[i];  
  }
}

filteredScalarField::filteredScalarField( std::string name, double* value)
:
 fieldName_(name),
 type(0)
{

 reference_=value;
 
 filter_size_[0]=-1;
 filter_size_[1]=-1;
 filter_size_[2]=-1;

}

filteredScalarField::filteredScalarField(std::string name,int number)
:
 fieldName_(name), 
 type(1)
{ 

 reference_=new double[number];
 
 filter_size_[0]=-1;
 filter_size_[1]=-1;
 filter_size_[2]=-1;

}

filteredScalarField::~filteredScalarField() 
{
if (type==1) delete reference_;  
}

