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

#ifndef C3PO_BASE_H
#define C3PO_BASE_H

#include "c3po.h"
#include "c3po_base_interface_vector.h"

namespace C3PO_NS
{

class c3poBase
{
  friend class c3poBaseAccessible;

  public:

    c3poBase(c3po *ptr) :
      pscl_(*ptr),
      c3poBaseInterfaceVector_(ptr->c3poBaseInterfaceVector_),
      operationContainer_(ptr->operationContainer_),
      selectorContainer_(ptr->selectorContainer_),
      control_(ptr->control_),
      comm_(ptr->comm_),
      input_(ptr->input_),
      output_(ptr->output_),
      error_(ptr->error_),
      dataStorageInternal_(ptr->dataStorageInternal_),
      timer_(ptr->timer_),
      mesh_(ptr->mesh_),
      basicMultiphaseQty_(ptr->basicMultiphaseQty_)
    {}

    virtual ~c3poBase() {}

  protected:

    inline c3po*                          c3po_ptr() const {return &pscl_;}
    inline const c3poBaseInterfaceVector& c3poInterface() const {return *c3poBaseInterfaceVector_;}

    inline const OperationContainer&      operationContainer() const {return *operationContainer_;}
    inline const SelectorContainer&       selectorContainer() const {return *selectorContainer_;}

    inline const Timer&                 timer()          const {return *timer_;}
    inline const Control&               control()        const {return *control_;}
    inline const Comm&                  comm()           const {return *comm_;}
    inline const Input&                 input()          const {return *input_;}
    inline const Output&                output()         const {return *output_;}
    inline const Error&                 error()          const {return *error_;}
    inline const DataStorage&           dataStorage()    const {return *dataStorageInternal_;}
    inline const c3poMesh&              mesh()                          const {return *mesh_;}
    inline const multiphaseFlowBasic&   basicMultiphaseQty()            const {return *basicMultiphaseQty_;}

  private:

    c3po &pscl_;
    c3poBaseInterfaceVector *&c3poBaseInterfaceVector_; // vehicle for global operations
    OperationContainer      *&operationContainer_;      // container of operations on data
    SelectorContainer       *&selectorContainer_;      // container of selectors
    Control                 *&control_;                 // control for stand-alone usage
    Comm                    *&comm_;                    // inter-processor communication
    Input                   *&input_;                   // input to simulation - input script processing
    Output                  *&output_;                   // input to simulation - input script processing
    Error                   *&error_;                    // error handling   
    DataStorage             *&dataStorageInternal_;     //data handling                                 
    Timer                   *&timer_;                   //timer performance recorder
    c3poMesh                *&mesh_;                    //mesh handling
    multiphaseFlowBasic     *&basicMultiphaseQty_;      //multiphase reference quantities
};

} //end c3po_NS

#endif

