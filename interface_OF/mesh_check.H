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

#ifndef mesh_check_C3PO_H
#define mesh_check_C3PO_H
#include "dictionary.H"
#include "coordinateSystem.H"
#include "coordinateSystems.H"
#include "wordList.H"
#include "labelList.H"
#include "dimensionedScalar.H"
#include "dimensionedTensor.H"
#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"
#include "fvMatricesFwd.H"
#include "fvMesh.H"
#include "c3po.h"
#include "error.h"
#include <mpi.h>
#include <vector>

namespace C3PO_NS
{

using namespace Foam;

class meshCheck
{
   private:
   
    const       Foam::fvMesh& mesh_;
    c3po*       C3po_;
    
     mutable int* centerID;
     mutable double mesh_length_[3]; 
     mutable double mesh_length_global_[3]; 	
     mutable double cell_middle_cell_0_[3];		
     mutable double cell_size_[3];
     mutable int    number_of_cells_[3];
     mutable int    number_of_cells_global_[3];

     mutable  double         minDomain[3];  
     mutable  double         maxDomain[3];
     mutable  double         minDomainGlobal[3];  
     mutable  double         maxDomainGlobal[3];
   public:
     meshCheck(const Foam::fvMesh&, c3po*);
     ~meshCheck();
     
     void checkMyMesh()  const;
     void determine_cellsize() const;
     void compute_globals() const;
     void create_centerID_vector()  const;
     void readBC() const;
     void registerC3POcells() const;

     int Cell_Number_x() {return number_of_cells_[0]; };
     int Cell_Number_y() {return number_of_cells_[1]; };
     int Cell_Number_z() {return number_of_cells_[2]; };
	
     int Cell_Number_global_x() {return number_of_cells_global_[0]; };
     int Cell_Number_global_y() {return number_of_cells_global_[1]; };
     int Cell_Number_global_z() {return number_of_cells_global_[2]; };
     
     double* min_Domain(){return minDomain; };
     double* max_Domain(){return maxDomain; };
     
     double* min_Domain_Global(){return minDomainGlobal; };
     double* max_Domain_GLobal(){return maxDomainGlobal; };
     
     double* cellSize() {return cell_size_;};


};

}
#endif
