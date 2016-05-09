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

#ifndef C3PO_BASE_ACCESSIBLE_H
#define C3PO_BASE_ACCESSIBLE_H

#include "c3po_base.h"

namespace C3PO_NS
{

class c3poBaseAccessible : c3poBase
{
  public:

    c3poBaseAccessible(c3po *ptr) : c3poBase(ptr) {}

    virtual ~c3poBaseAccessible() {}

  protected:

    inline c3po*        c3po_ptr()    const {return &pscl_;}
    inline OperationContainer& operationContainer() {return *operationContainer_;}
    inline SelectorContainer&  selectorContainer() {return *selectorContainer_;}
    inline DataStorage&    dataStorage() const {return *dataStorageInternal_;}
    inline Control&        control()  const {return *control_;}
    inline Comm&           comm()     const {return *comm_;}
    inline Input&          input()    const {return *input_;}
    inline Output&         output()   const {return *output_;}
    inline Error&          error()    const {return *error_;}
    inline Timer&          timer()    const {return *timer_;}
    inline c3poMesh&       mesh()     const {return *mesh_;}
    inline multiphaseFlowBasic&  basicMultiphaseQty()     const {return *basicMultiphaseQty_;}

  private:

};

} //end c3po_NS

#endif

