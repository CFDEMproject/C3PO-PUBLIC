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
	Base class for a selector (i.e., a class that manages the selection of cell
    or particle indices).
-----------------------------------------------------------------------------------*/

#ifndef C3PO_SELECTOR_BASE_H
#define C3PO_SELECTOR_BASE_H

#include "c3po_base_accessible.h"
#include "string.h"
#include "input.h"
#include "error.h"
#include "data_storage.h"
#include "mesh.h"
namespace C3PO_NS
{

//numbering of Selector type
enum{ CELL,		//0
      PARTICLE,   	//1
      REGION            //2
    };


class SelectorBase : public c3poBaseAccessible
{
    public:

      SelectorBase(c3po *ptr,const char *_name);
      ~SelectorBase();

      virtual void init(int narg, char const* const* arg) {};

      virtual void begin_of_step() {};
      virtual void middle_of_step() {};
      virtual void end_of_step() {};

      virtual int findNearestCell(double x, double y, double z) {return -1;};
      
      //Access functions
      const char* name() {return name_;}
      const char* scope() {  return scope_; }
      virtual double value(){return 0.0;};
      
      
      inline double* currentCellCenters() {return cellList_;}; 
      
      bool haveCell() {
                       if( currentCell_<mesh().MaxNofCellsProc()[comm().me()]
                           &&
                           inside_[comm().me()]
                         )
                        return true;
                        else return false;
                       }
      bool haveParticle() {
                       if( currentPar_<dataStorage().getNofParticlesProc(comm().me()) 
                           &&
                           inside_[comm().me()]
                         )
                        return true;
                        else return false;
                       }
      

    private:

      char *name_;
      char *scope_; // scope/name for JSON file
    
    protected:
      
       int currentCell_;
       int currentPar_;
       bool* inside_;
       double* cellList_;
};

} //end c3po_NS

#endif

