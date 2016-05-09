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
#include "lagrangian.h"
#include "c3po.h"
#include <string>
#include <iostream>
#include <stdlib.h> 
#include <cmath>
#include  <string>
#include  <fstream>
#include  <vector>
#include <cstdlib>

#define  MAX_CHARS_PER_LINE3 1024

using namespace C3PO_NS;

/*-----------------------------------------------------------------------------*/
CSVlagrangian::CSVlagrangian(c3po * c3po_, CSVmesh * mesh)
:
C3po_(c3po_),
mesh_(mesh),
NofPar_(0)
{
 tolerance_ = C3po_->meshCheckTolerance();
} 

/*-----------------------------------------------------------------------------*/
CSVlagrangian::~CSVlagrangian()
{
 deleteParticles();
}
/*-----------------------------------------------------------------------------*/
void CSVlagrangian::checkLagrangian() const
{
 std::ifstream infile_ ("lagrangian.csv",std::ifstream::in);
    
 if (!infile_.is_open())
 {
  std::cout << "\n-> WARNING: No lagrangian.csv file was found! Particles will not be registered\n";
  return; 
 }

 readLagrangian();
 registerC3POparticles(); 
 


}
/*-----------------------------------------------------------------------------*/
void CSVlagrangian::readLagrangian() const
{
  std::ifstream infile_ ("lagrangian.csv",std::ifstream::in);
  std::istream& infile(infile_);  
  
 int s=0;
 std::vector<int> map;
 
 int decDir_ = C3po_->getDomainDecomposition();

 while (!infile.eof())
 {
  
  //read the line
  std::vector<double> temp;
  std::string buf;
  std::getline(infile,buf);
  std::string::iterator it=buf.begin();
  while(it<buf.end())
  {
   std::string str;
   while(it!=buf.end())
   {
    if(*it==',') break;
    str.push_back(*it);
    it++;
   }
     
   temp.push_back(std::atof(str.c_str()));  
   it++;
  }
  
  if(temp.size()==0) break;
  //check if inside the local domain
  if(   ( temp[decDir_] + tolerance_ ) > mesh_->localMin()[decDir_] 
     && ( temp[decDir_] - tolerance_ ) < mesh_->localMin()[decDir_]
    )
  {
   
   int size_ = temp.size();
   
   if (size_ < 13 )
   {
    std::cout << "\n -> ERROR: You have to enter, at least: radius, position (x,y,z), velocity (x,y,z), torque (x,y,z) and on force (x,y,z) in lagrangian.csv (13 columns)\n";
    MPI_Finalize();
    exit(1);
   }
   
   particleRadius_.push_back(temp[0]);
   
   double * tempcoord_   = new double[3];
   double * tempvel_     = new double[3];
   double * tempTorque_ = new double[3];
   for(int i=0;i<3;i++)
   {
    tempcoord_[i] = temp[i+1];
    tempvel_[i]   = temp[i+4];
    tempTorque_[i] = temp[i+7];
   }
   particlePos_.push_back(tempcoord_);
   particleVel_.push_back(tempvel_);
   particleTorque_.push_back(tempTorque_);
   NofPar_++;
   
   int forces_ = (size_ - 10) / 3; 
   
   std::vector<double*> tempVec_;
   
   for(int i=0;i<forces_;i++)
   {
    double * tmpforce_ =new double[3];
    
    tmpforce_[0]= temp[10 + 3*i];
    tmpforce_[1]= temp[11 + 3*i];
    tmpforce_[2]= temp[12 + 3*i];
    
    tempVec_.push_back(tmpforce_);  
    
    NofPar_++;
    
   }
   
   particleForce_.push_back(tempVec_);
   
  }
  s++;
 } 
 

}
/*------------------------------------------------------------------------------*/
void CSVlagrangian::registerC3POparticles() const
{
 for(int i=0;i<NofPar_;i++)
   C3po_->addParticle("CSVparticles", particleRadius_[i], particlePos_[i], particleVel_[i], &particleForce_[i],NULL, particleTorque_[i]);

}
/*------------------------------------------------------------------------------*/
void CSVlagrangian::deleteParticles() const
{
 for(int i=0;i<NofPar_;i++)
 {
  delete particlePos_[i];
  delete particleVel_[i];
  delete particleTorque_[i];
  
  for(unsigned int j=0;j<particleForce_[i].size();j++)
   delete particleForce_[i][j];
     
 }
 particlePos_.clear(); 
 particleVel_.clear(); 
 particleForce_.clear(); 
 particleTorque_.clear(); 
 particleRadius_.clear(); 
  
 NofPar_=0;
}
