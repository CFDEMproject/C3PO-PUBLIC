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
#include "mesh_check.h"
#include "c3po.h"
#include <string>
#include <iostream>
#include <stdlib.h> 
#include <cmath>
#include  <string>
#include  <fstream>
#include  <vector>
#include <cstdlib>

#define  MAX_CHARS_PER_LINE2 1024

using namespace C3PO_NS;

CSVmesh::CSVmesh(c3po * c3po_)
:
C3po_(c3po_),
NofCells_(0),
nprocs_(1),
me_(0)
{
 MPI_Comm_rank(MPI_COMM_WORLD,&me_);
 MPI_Comm_size(MPI_COMM_WORLD,&nprocs_);
 tolerance_ = C3po_->meshCheckTolerance();
}
/*-----------------------------------------------------------------------------*/
CSVmesh::~CSVmesh()
{
 clearMesh();
}
/*-----------------------------------------------------------------------------*/
void CSVmesh::checkMesh() const
{
 computeParallelDomain();
 readMeshLocal();
 registerC3poCells();
 
}
/*-----------------------------------------------------------------------------*/
void CSVmesh::computeParallelDomain() const
{

 for(int i=0;i<3;i++)
 {
  MaxGlobalDomain_[i] = C3po_->getMaxDomainfromJson(i);
  MinGlobalDomain_[i] = C3po_->getMinDomainfromJson(i);
 }

 decDir_ = C3po_->getDomainDecomposition();
 
 double delta_;
 

 delta_ = ( MaxGlobalDomain_[decDir_] - MinGlobalDomain_[decDir_] ) / nprocs_;
  
 for(int i=0;i<3;i++)
 {
  if(i == decDir_ && nprocs_ > 1) continue;
  MaxLocalDomain_[i]=MaxGlobalDomain_[i];
  MinLocalDomain_[i]=MinGlobalDomain_[i];
 }
 
 if(nprocs_ > 1)
 {
  MinLocalDomain_[decDir_] = me_ * delta_ + MinGlobalDomain_[decDir_];
  MaxLocalDomain_[decDir_] = MinLocalDomain_[decDir_] + delta_;
 } 
 //register in CPPPO
 
 
}

/*-----------------------------------------------------------------------------*/
void CSVmesh::readMeshLocal() const
{
 std::ifstream infile_ ("mesh.csv",std::ifstream::in);
 int columnsToRead[4];
 
 C3po_->csv_columnsToRead(&columnsToRead[0]);
    
 if (!infile_.is_open())
 {
  std::cout << "\n proc: " << me_ << " -> ERROR: mesh.csv can not be opened!\n";
  MPI_Finalize();
  exit(1);
 }

 std::istream& infile(infile_);    
 
 //Number of processed cells 
 int s=0;
 
 //Get rid of the header
 std::string header;
 std::getline(infile,header);
 
 while (!infile.eof())
 {
    //read the line
  double temp[4];
  std::string buf;
  std::getline(infile,buf);
  std::string::iterator it=buf.begin();
  int colcount=0;
  int column = 0;
  
  while(it<buf.end())
  {
  
   
   std::string str;
   
   while(it!=buf.end())
   {
    if(*it==',') break;
    if(*it==' ')
    {
     it++;
     continue;
    }
    str.push_back(*it);
    it++;
   }
   
   for(int var=0;var<4;var++)  
   {
    if(column ==  columnsToRead[var])
    { 
     temp[var]= std::atof(str.c_str());  
     colcount++;
    }
   }
   
   it++;
   column++;
  }
  
  if(colcount==0 && s>0) break;
  if(colcount != 4) 
  {
   std::cout << "\n proc: " << me_ << " -> ERROR: problems with the mesh file! Be sure to provide the four columns in the correct format!\n";
   MPI_Finalize();
   exit(1);
  } 
  
  //check if inside the local domain
  if(   ( temp[decDir_] + tolerance_ ) > MinLocalDomain_[decDir_] 
     && ( temp[decDir_] - tolerance_ ) < MaxLocalDomain_[decDir_]
    )
  {
   
   for(int i=0;i<3;i++)
    cellCoord_.push_back(temp[i]);
   
   cellVol_.push_back(temp[3]);
   NofCells_++;
   cellId_.push_back(s);
  }
  s++;
 } 
 
 NofCellsGlobal_=s;
  //std::cout << "\n proc: " << me_ << " -> Number of my cells found in mesh.csv: "<< NofCells_ << " \n";
 
}
/*-----------------------------------------------------------------------------------*/
void CSVmesh::registerC3poCells() const
{
 C3po_->setNofCells(NofCells_);
 
 C3po_->registerCells(&cellVol_[0] ,&cellCoord_[0],NofCells_,3) ;
  
 C3po_->registerDomainInfo(MaxLocalDomain_,MinLocalDomain_,MaxGlobalDomain_,MinGlobalDomain_);
}
/*-----------------------------------------------------------------------------------*/
void CSVmesh::clearMesh() const
{
 
 cellCoord_.clear();
 cellVol_.clear();
 cellId_.clear();
 
 NofCells_ = 0;
}

