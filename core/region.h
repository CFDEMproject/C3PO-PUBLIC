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
     Class for holding region information (e.g., to be used in a selector
-----------------------------------------------------------------------------------*/

#ifndef C3PO_region_H
#define C3PO_region_H


#include "c3po_base_accessible.h"
#include "comm.h"
#include <vector>
#include "qjson_includes.h"
#include <string>

namespace C3PO_NS
{
 class region : public c3poBaseAccessible
 {
  private:
  
  std::string                             regionName_;    //identifies the region
  
  int                                       regionId_;      //identifies the region

  std::vector<std::string>               filterNames_;   //name of the filters that use this region
  
  mutable std::vector<int>           cellIDsInternal_; //Holds ids of internal cells (Note that also surface cells are included!)

  mutable std::vector<int>            cellIDsSurface_;  //Holds ids of surface cells (These are just surface cells)
 
  mutable int*                        nParticlesProc_;
  
  mutable bool                     readCellsFromJson_;
  
  double                                      volume_;
  
  double                                 surfaceArea_;
  
  double                             centerOfMass_[3];
  
  public:
 
 
  region(std::string regionName, std::vector<std::string> filterNames, bool readCellsFromJson, c3po *ptr);
  ~region();

  //Functionality 
  void addCell(int & IDInternal, bool isOnSurface);

  void clearCells();

  void evaluateVolume();

  void evaluateSurfaceArea();

  void evaluateCenterOfMass();
  
  void write();
  
 
  //Access  
  int numOfCells() {return cellIDsInternal_.size();};

  int numOfCellsSurface() {return cellIDsSurface_.size();};

  double volume() { return volume_;};

  double surfaceArea() {return surfaceArea_;};
  
  double centerOfMass(int i) {return centerOfMass_[i];};
   
 };
}



#endif
