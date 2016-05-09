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
    Class to parse and store mesh data from CSV.
-----------------------------------------------------------------------------------*/
#ifndef CSV_MESH_CHECK
#define CSV_MESH_CHECK
#include "c3po.h"
#include  <string>
#include  <fstream>
#include  <vector>

namespace C3PO_NS
{

 class CSVmesh
 {
  
  public:
  
  CSVmesh(c3po*);
  ~CSVmesh();
  
  void checkMesh() const;
   
  inline int* getCellId(int i) const {return &(cellId_[i]);};
  inline int* NofCells() const {return &NofCells_;}; 
  inline int* NofCellsGlobal() const {return &NofCellsGlobal_;}; 
  
  inline double* localMax() const {return &MaxLocalDomain_[0];};
  inline double* localMin() const {return &MinLocalDomain_[0];};
  
  void clearMesh() const;
  
  
  private:
  
  c3po* C3po_;
  
  void computeParallelDomain() const;
  void readMeshLocal() const;
  
  void registerC3poCells() const;
  
  
  
  
  mutable double tolerance_;
  
  mutable std::vector<double>  cellCoord_;
  mutable double               MaxLocalDomain_[3];
  mutable double               MinLocalDomain_[3];
  
  mutable double               MaxGlobalDomain_[3];
  mutable double               MinGlobalDomain_[3];
  
  mutable int NofCells_;
  mutable int NofCellsGlobal_;
  
  mutable std::vector< double > cellVol_;
  
  mutable int nprocs_;
  mutable int me_;
  
  mutable int decDir_;
  
  mutable std::vector<int> cellId_;
 
 };
}


#endif
