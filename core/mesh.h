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
	Class to hold references and functions to manage mesh info. 
-----------------------------------------------------------------------------------*/

#ifndef C3PO_MESH_H
#define C3PO_MESH_H

#include "stdio.h"
#include "c3po_base.h"
#include "c3po_base_interface.h"
#include "output.h"
#include "comm.h"
#include <map>
#include <string>
#include "qjson_includes.h"

#include "data_storage.h"


namespace C3PO_NS
{

class c3poMesh : public c3poBase, public c3poBaseInterface
{
  public:

  c3poMesh(c3po *ptr);
  ~c3poMesh();
      
     
 double meshCheckTolerance() const;
 double meshFilterWidthTolerance() const;
 bool   meshVerbose() const;
      
 void    setCellsize(double* x) const  {
	     cellSize_[0]=x[0]; 
	     cellSize_[1]=x[1];
	     cellSize_[2]=x[2]; 
	     haveCellSize_=true;
	  };
     
 double *cellSize() const {  return cellSize_;};
     
 bool    haveCellSize() const {return haveCellSize_;};
   
 inline int* cell_global_x()	    const {return &CellCount_global_[0];};  
 inline int* cell_global_y()	    const {return &CellCount_global_[1];};	
 inline int* cell_global_z()         const {return &CellCount_global_[2];};
	
 inline int* cell_global()           const {return &CellCount_global_[0];}; 
	
 inline int totalGlobalCell()        const {return CellCount_global_[0]*CellCount_global_[1]*CellCount_global_[2];};	
	
 inline int* cell_x()		    const {return &CellCount_[0];};  
 inline int* cell_y()		    const {return &CellCount_[1];};	
 inline int* cell_z()		    const {return &CellCount_[2];};
	
 inline int* cells_subdomain()       const {return &CellCount_[0];}; 
	
 inline int  NofCells()              const {return NofCells_;};	
 	          
 void gatherInfo() const;
	
 //Functions for managing unstructured mesh data
// void registerCell(double cellV,const double* centreX_,const double* centreY_,const double* centreZ_) const;  
 void clearCells() const;
	
 inline double* CellCentre(int id)     const {return  &cellCentre_[meshXYZspacing_*id];}; 
	
 inline double*  CellVol(int id)        const {return &cellV_[id];};
	
 inline double* dMax()                const {return complete_domain_max;};
 inline double* dMin()                const {return complete_domain_min;};
	
 inline double* dLocalMax()                const {return &local_domain_max[0];};
 inline double* dLocalMin()                const {return &local_domain_min[0];};
	
 inline double* dGlobalMax()           const {return &global_domain_max[0];};
 inline double* dGlobalMin()           const {return &global_domain_min[0];};
		
 void  registerDomainInfo(double maxDomain[3],double minDomain[3],double maxDomainGlobal[3],double minDomainGlobal[3]) const;    
	
 void setNofCells(const int x) const {NofCells_=x;};
 
 void setcell(int* CellCount) const;
 
 void setcell_global(int *global_number_of_cells) const;
	
 inline int* MaxNofCells()     const {return &MaxNofCells_;};
 inline int* MaxNofCellsProc()     const {return MaxNofCellsProc_;};
	
 inline int* ijkMax()              const {return  complete_IJKdomain_max;};
 inline int* ijkMin()              const {return  complete_IJKdomain_min;};
	
 inline int* localIJKMax()              const {return  &local_ijk_max[0];};
 inline int* localIJKMin()              const {return  &local_ijk_min[0];};
 
 void calculateIJKfromJson() const;
 
 void setMeshXYZSpacing(int space_)              {meshXYZspacing_=space_;};   
 
 void registerCells(double* Vvector, double* posVector, int NumberOfCells, int meshXYZspacing);
	                   
 private:      
 
 double*                                      cellV_;
 double*                                 cellCentre_;
 
 int                                 meshXYZspacing_;
           
 mutable double*   complete_domain_max;
 mutable double*   complete_domain_min;
      
 mutable double    global_domain_max[3];
 mutable double    global_domain_min[3];
      
 mutable double    local_domain_max[3];
 mutable double    local_domain_min[3];

 mutable int*      complete_IJKdomain_max;
 mutable int*      complete_IJKdomain_min;
      
 mutable int       local_ijk_max[3];
 mutable int       local_ijk_min[3];
            
 mutable int NofCells_;
 mutable int MaxNofCells_;
 mutable int* MaxNofCellsProc_;
      
 mutable int CellCount_global_[3]; 
 mutable int CellCount_[3];                	

 mutable double cellSize_[3];  //note, all cells must have the same size
 mutable bool   haveCellSize_;
	  

};

} //end c3po_NS

#endif
