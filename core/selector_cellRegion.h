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
	Selector class for region detection and handling
-----------------------------------------------------------------------------------*/

#ifdef SELECTOR_CLASS

SelectorStyle(bubble, SelectorCellRegion)

#else


#ifndef C3PO_SELECTOR_CELLREGION_H
#define C3PO_SELECTOR_CELLREGION_H

#include "selector_base.h"

namespace C3PO_NS
{


class SelectorCellRegion : public SelectorBase
{
    public:

      SelectorCellRegion(c3po *ptr,const char *_name);
      ~SelectorCellRegion();

//      void init(int narg, char const* const* arg) {};

      void begin_of_step();
      void middle_of_step() ;
      void end_of_step() {};

    private:

      //Variable
      double                threshold_;                 //for detecting the region
      bool                  identifyBelowThreshold_;    //for detecting the region
      std::string           nameOfFieldForSelection_;   
      int                   fieldForSelection_;         //pointer to field
      bool                  updateRegionsEachStep_;     //if true, every time step the regions will be updated
      bool                  updated_;
      int*                  AA;                        //bubble vector. It tells if a bubble is present.
      int                   NofBubbles_;
 
      //Functionality
      void (SelectorCellRegion::*updateRegions)();
      
      void IJKRegionSelector();
      
      int ijk2CellIDOF(int i_add_correct, int j_add_correct, int k_add_correct);
      void CellID2ijk(int center_Cell_ID, int* result); 
      
      bool (SelectorCellRegion::*check)(int);
      
      bool checkIfBelow(int cell);
      bool checkIfAbove(int cell);
      
      std::vector<int>  bubble_index_store;
      std::vector<int> bubble_index_surface;
      int               bubble_index_eval;
      
      void checkIndexStore(int id);
      void createRegion();
    
};

} //end c3po_NS

#endif
#endif
