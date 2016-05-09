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
/*-----------------------------------------------------------------------------------
Description
	This class holds the references and the functions needed to access data from
	a single particle.
-----------------------------------------------------------------------------------*/

#ifndef PARTICLE_H
#define PARTICLE_H
#include<vector>

namespace C3PO_NS
{
 class Particle
 {
  private:
  double  buf[3];
  
  int id_;
  
  double                  r_;
  double*                 pos_;
  double*                 vel_;
  double*                 torque_;
  std::vector< double* >* force_;
  std::vector< double >*  scalar_;
  
  int                     cellCentreId_;
  
  std::vector<double*>    filteredVectors_;
  std::vector<double>     filteredScalars_;
  std::vector<double* >   forceBuf_[1];
  std::vector<double>     scalarBuf_;
  
  double                  filterVolume_;
  
  
  public:
  
  Particle();
  ~Particle();
  
  void setradius(double x)                 {r_=x;};
  void setpos(double* x )                  {pos_=x;};
  void setvel(double* x)                   {vel_=x;};
  void settorque(double * x)               {torque_=x;};
  void setforce(std::vector< double* >* x) {force_=x;};
  void setscalar(std::vector< double >* x) {scalar_=x;};
  void setCellCentreId_(int x)             {cellCentreId_=x;};    
  void setId(int x)                        {id_=x;};
  void setFilterVolume(double x)            {filterVolume_=x;};
  
  inline double*    getradius()            {return &r_;};
  inline double* getpos()                  {return pos_;};
  inline double* getvel()                  {return vel_;};
  inline double* gettorque()               {return torque_;};
  inline double* getforce(int i)           {return (*force_)[i];};
  inline double  getscalar(int i)          {return (*scalar_)[i];};
  inline int getNofForces()                {return force_->size();};
  inline int getId()                       {return id_;};
  inline double getFilterVolume()          {return filterVolume_;};
  
  int   checkScalars();                    
  
  double getTotalForce(int i);
  
  int returnCellCentreId_()               {return cellCentreId_;};

  inline double*  filteredVector(int i)  {return filteredVectors_[i];};
  inline double*  filteredScalar(int i)   {return &(filteredScalars_[i]);};
  
  int NofVF() {return filteredVectors_.size();};
  int NofSF() {return filteredScalars_.size();};
  
  void  addVector();
  void  addScalar(); 
  
  void deleteContent();
 };
}
#endif
