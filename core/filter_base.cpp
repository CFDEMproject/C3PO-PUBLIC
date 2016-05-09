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

#include "selector_container.h"
#include "filter_base.h"
#include <string>
#include "comm.h"
#include "mesh.h"
#include "output.h"
#include "error.h"
#include <fstream>
#include <cmath>
#define PI 3.14159265359

using namespace C3PO_NS;
using std::ifstream;

/* ----------------------------------------------------------------------
   FilterBase Constructor
------------------------------------------------------------------------- */

FilterBase::FilterBase(c3po *ptr,const char *name, int myType) 
: 
    c3poBaseAccessible(ptr),
    typeRegion_(myType),
    name_(0),
    cellSize_(NULL),
    total_filtered_cells_(0),
    selective_(false)
{
    int n = strlen(name) + 1;
    name_ = new char[n];
    strcpy(name_,name);
    
    
    cellSize_ = mesh().cellSize();
}

FilterBase::~FilterBase()
{
    delete name_;
    
}


/* ----------------------------------------------------------------------
   Init with QJson data
------------------------------------------------------------------------- */
void FilterBase::init(QJsonObject jsonObj) const
{

  
    if (jsonObj["CoordSys"].isNull())
        error().throw_error_one(FLERR,"You must specify a coordinate system for each filter. \n");

    typeRegion_    = jsonObj["CoordSys"].toInt();
    if (jsonObj["CoordSys"].toInt()==0)
    {
     filtersize_[0] = (jsonObj["x"].toDouble()-c3po_ptr()->meshFilterWidthTolerance())/2; //avoid problems in case filter width matches cell-cell dist.
     filtersize_[1] = (jsonObj["y"].toDouble()-c3po_ptr()->meshFilterWidthTolerance())/2;
     filtersize_[2] = (jsonObj["z"].toDouble()-c3po_ptr()->meshFilterWidthTolerance())/2;
    }
    else if(jsonObj["CoordSys"].toInt()==1)
    {
     filtersize_[0]= jsonObj["r"].toDouble()+c3po_ptr()->meshFilterWidthTolerance();
     filtersize_[1]= jsonObj["r"].toDouble()+c3po_ptr()->meshFilterWidthTolerance();
     filtersize_[2]= jsonObj["r"].toDouble()+c3po_ptr()->meshFilterWidthTolerance();
    }

    filterType_ = jsonObj["CoordSys"].toInt();
    
    if(!jsonObj["selective"].isNull())
    {
     
     if(jsonObj["selective"].toBool())
     {
      selective_=true;
      if(jsonObj["max"].isNull())
       error().throw_error_one(FLERR,"You must specify the max domain size when using selective filters. \n"); 
      
      if(jsonObj["min"].isNull())
       error().throw_error_one(FLERR,"You must specify the min domain size when using selective filters. \n"); 
      
      for(int i=0;i<3;i++)
      {
       maximum_[i]=jsonObj["max"].toArray()[i].toDouble();
       minimum_[i]=jsonObj["min"].toArray()[i].toDouble();
      }
      
     }    
    }
  
}

void FilterBase::setNFilterCells() const
{
    if(typeRegion_==CARTESIAN)
    {
    
    //Only proc 0 does this calculation and, then, broadcast...this to prevent processors from reading different filter sizes.
     if(comm().is_proc_0())
     {
      for(int iDir=0; iDir<3; iDir++)
       filterNcells_[iDir] = lround((filtersize_[iDir] )/(cellSize_[iDir]+ c3po_ptr()->meshFilterWidthTolerance()));
     } 
    MPI_Barrier(MPI_COMM_WORLD);
   
    for(int iDir=0; iDir<3; iDir++)
     MPI_Bcast(&filterNcells_[iDir],sizeof(int),MPI_INT,0,MPI_COMM_WORLD);

    if(input().verbose())
     {
      total_filtered_cells_ = (filterNcells_[0]*2+1)*(filterNcells_[1]*2+1)*(filterNcells_[2]*2+1);
      printf("type: %d, filtersize = %gx%gx%g (%dx%dx%d cells) on grid with cellSize %gx%gx%g results in %d filtered cells. \n", 
                 typeRegion_,
                 filtersize_[0], 
                 filtersize_[1],
                 filtersize_[2],
                 filterNcells_[0],
                 filterNcells_[1],
                 filterNcells_[2],
                 cellSize_[0],
                 cellSize_[1],
                 cellSize_[2],
                 total_filtered_cells_
                 
    
    
          );
      }
    }
    else if(typeRegion_==SPHERICAL)
    {
    
   //Only proc 0 does this calculation and, then, broadcast...this to prevent processors from reading different filter sizes.
     if(comm().is_proc_0())
     {
      for(int iDir=0; iDir<3; iDir++)
       filterNcells_[iDir] = lround((filtersize_[iDir] -cellSize_[iDir]/2)/(cellSize_[iDir]+ c3po_ptr()->meshFilterWidthTolerance()));
     } 
    MPI_Barrier(MPI_COMM_WORLD);
   
    for(int iDir=0; iDir<3; iDir++)
     MPI_Bcast(&filterNcells_[iDir],sizeof(int),MPI_INT,0,MPI_COMM_WORLD);

    if(input().verbose())
     {
      total_filtered_cells_ = (filterNcells_[0]*2+1)*(filterNcells_[1]*2+1)*(filterNcells_[2]*2+1);
      printf("type: %d, filtersize = %gx%gx%g (%dx%dx%d cells) on grid with cellSize %gx%gx%g results in %d filtered cells. \n", 
                 typeRegion_,
                 filtersize_[0], 
                 filtersize_[1],
                 filtersize_[2],
                 filterNcells_[0],
                 filterNcells_[1],
                 filterNcells_[2],
                 cellSize_[0],
                 cellSize_[1],
                 cellSize_[2],
                 total_filtered_cells_
                 
    
    
          );
      }
    }       
    else
        error().throw_error_one(FLERR,"This region type is not implemented! \n");
 if(input().verbose())
    printf("nFilter with name %s initialized and filter cell count set. \n\n", name_);

}
