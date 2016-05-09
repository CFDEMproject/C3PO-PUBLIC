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

#include "particle.h"
#include <cstddef>
using namespace C3PO_NS;

Particle::Particle()
:
id_(-1)
{
 buf[0]=0.0;
 buf[1]=0.0;
 buf[2]=0.0;
 torque_=&buf[0];
 vel_=&buf[0];
 pos_=&buf[0];
 cellCentreId_=-1;
 force_=&forceBuf_[0];
 scalar_=&scalarBuf_;
}

Particle::~Particle() 
{
 deleteContent();
}
/* * * * * * * * * * * * * * * * * * * * * * * * * * */
double Particle::getTotalForce(int i)
{
 double totalForce_=0;
 for(unsigned int j=0; j<force_->size(); j++  )
  totalForce_+=(*force_)[j][i];
 
 return totalForce_;

} 
/* * * * * * * * * * * * * * * * * * * * * * * * * */
void Particle::addVector()
{
  double * tmp_ = new double[3];
  
  for(int i=0;i<3;i++)
   tmp_[i]=.0;
   
  filteredVectors_.push_back(tmp_);
  
}
/* * * * * * * * * * * * * * * * * * * * * * * * * */
void Particle::addScalar()
{
   
  filteredScalars_.push_back(0);
  
}
/* * * * * * * * * * * * * * * * * * * * * * * * * */
void Particle::deleteContent()
{
  for(unsigned int i=0;i<filteredVectors_.size();i++)
   delete filteredVectors_[i];
   
  filteredVectors_.clear();
  filteredScalars_.clear();
  
}
/* * * * * * * * * * * * * * * * * * * * * * * * * */
int Particle::checkScalars()
{
 return (scalar_== NULL )? -1 : scalar_->size();
};

