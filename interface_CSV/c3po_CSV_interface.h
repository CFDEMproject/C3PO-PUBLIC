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
    Main class for CPPPO-CSV interface.
-----------------------------------------------------------------------------------*/

#ifndef c3poCSVInterface_H
#define c3poCSVInterface_H

#include  "c3po.h"
#include  <string>
#include  <fstream>
#include  <vector>
#include "mesh_check.h"
#include "lagrangian.h"
#include "CSVfieldOperations.h"

namespace C3PO_NS
{
 class c3poCSVinterface
 {
  public:
  
  c3poCSVinterface(MPI_Comm);
  ~c3poCSVinterface();
  
  void checkMesh() const;
  
  void clearMesh() const;
  
  void checkParticles() const;
  
  void clearParticles() const;
    
  void runC3po() const;
  
  void FLUENTrunC3PO() const;
    
  private:
  
  void parseFile() const;
  
  void registerC3poFields() const; 
  
  void printCSV(int id) const;
  
  void readInput() const;
  
  void createFileList() const;
  
  void deleteC3POfields() const;
  
  void runFilter(int id) const;
  
  void runSampling() const;
  
  void runBinning() const;
  
  void createFilterFields(int id) const;
  
  void deleteFields() const;
  
  void resetAllFields() const;
  
  void createGradients(int id) const;

  int check_repetition (int n, int sf) const; 
  
  
  c3po* myC3po_;
  
  CSVmesh * mesh_;
  
  CSVlagrangian * lagrangian_;
  
  mutable std::string fileName_;
  
  MPI_Comm comm_;
  
  int comm_me_;
  
  int nprocs_;
  
  mutable int        NofCells_;
  
  mutable int       NofFields_;
  
  mutable int      NofVectors_;
  
  mutable double**        csv_;
  
  mutable double**     fields_;
  
  mutable std::vector<std::string>  names_;
  
  mutable std::vector<std::string>  Ffieldnames_;
  
  mutable char**      vecName_;
  
  mutable int   NofFiltFields_;
  
  mutable int   NofFiltVarianceFields_;
  
  mutable int   timeId_;
  
  mutable bool  outDirGenerated_;
  
  mutable bool twoD_;
  
  mutable std::vector<double*> dummyZ;
  
  mutable std::vector<std::string> fileList_;
  
  mutable std::vector<double*>          gradientFields_;
  
  CSVfieldOperations*              fO_;
  
 };
}

#endif
