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
	The base (default) class to specify the size of the filter region
-----------------------------------------------------------------------------------*/

#ifndef C3PO_FILTER_BASE_H
#define C3PO_FILTER_BASE_H

#include "c3po_base_accessible.h"
#include <QtCore/QJsonObject>

#define SMALLNUMBER 1e-32

namespace C3PO_NS
{

//numbering of filter region type
enum{ CARTESIAN, 	//0 - useful for filtering on eulerian grids, default
      SPHERICAL,        //1 - useful for filtering particle data, 
      CYLINDRICAL       //2 - useful for filtering particle data, not implemented
    };

class FilterBase : public c3poBaseAccessible
{
    public:

      FilterBase(c3po *ptr,const char *_name, int myType);
      ~FilterBase();

      virtual void init(QJsonObject jsonObj) const;

      void   setNFilterCells() const;
      
      //Access functions
      const char* name() {return name_;};
      int filterNcells(int x) {return filterNcells_[x];};
      int filterType() {return filterType_;};
      inline double* filterSize(int x)  {return &filtersize_[x];};
      
      bool selectiveFilter()          {return selective_;};
      
      double* maxToFilter()           {return &maximum_[0];};
      double* minToFilter()           {return &minimum_[0];};

    private:

      mutable int  typeRegion_;

      char *name_;

      //for CARTESIAN type
      double          *cellSize_;            //spatial extension of a single cell
      mutable double  filtersize_[3];        //spatial extension of filter
      mutable int     filterNcells_[3];      //cell count to filter over
      mutable int     total_filtered_cells_; //total number of cells to filter
      mutable int     filterType_;
      mutable bool    selective_;            //true if just a portion of the domain has to be filtered
      
      mutable double  minimum_[3];           //minimum domain to filter
      mutable double  maximum_[3];           //maximum domain to filter
      
};

} //end c3po_NS

#endif

