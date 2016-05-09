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
  This class contains the functions needed to select one or more cells (position to cellId)
  using an OpenFOAM general unstructured mesh. Only a linear finding algorithm is implemented. 
-----------------------------------------------------------------------------------*/

#ifdef SELECTOR_CLASS

SelectorStyle(cellUnstruct, SelectorCellUnstruct)

#else


#ifndef C3PO_SELECTOR_CELL_UNSTRUCT_H
#define C3PO_SELECTOR_CELL_UNSTRUCT_H

#include "c3po_base_accessible.h"
#include "selector_base.h"

namespace C3PO_NS
{


class SelectorCellUnstruct : public SelectorBase
{
    public:

      SelectorCellUnstruct(c3po *ptr,const char *_name);
      ~SelectorCellUnstruct();

      void begin_of_step();
      void middle_of_step() {};
      void end_of_step() {};
      
      int findNearestCell(double x, double y, double z);

      
    private:
    
    double tolerance_;
    
    

    double* max_;
    double* min_;
    bool periodic_[3];

    bool (SelectorCellUnstruct::*checkX)(double,double,double);
    bool (SelectorCellUnstruct::*checkY)(double,double,double);
    bool (SelectorCellUnstruct::*checkZ)(double,double,double);
    
    bool (SelectorCellUnstruct::*checkR)(double,double,double,double,int,double*);
    
    bool checkPosition(double posA, double max, double min );
    bool checkPositionPeriodic(double posA, double max, double min );
    
    bool checkPositionSphere(double i, double j, double k, double r_, int index_, double * delta_);
    bool checkPositionPeriodicSphere(double i, double j, double k, double r_, int index_, double * delta_);
    bool checkPositionSphereDummy(double i, double j, double k, double r_, int cell_, double * delta_) {return true;};
        
    void RunSelector();
    void fillArrays();
    void boundaryCorrections();
    
    
    

};

} //end c3po_NS

#endif
#endif
