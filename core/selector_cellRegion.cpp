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


#include "selector_cellRegion.h"
#include "input.h"
#include "mesh.h"
#include "error.h"
#include "qjson_includes.h"
#include "data_storage.h"

#define TOLERANCE 1E-10


using namespace C3PO_NS;

/* ----------------------------------------------------------------------
   ModelBase Constructor
------------------------------------------------------------------------- */

SelectorCellRegion::SelectorCellRegion(c3po *ptr,const char *name) 
: 
SelectorBase(ptr,name),
identifyBelowThreshold_(true),
updateRegionsEachStep_(false),
updated_(false),
NofBubbles_(0)
{
   if(input().mainDoc()->object()[name].isNull())
            {
             printf("\nERROR: can not find the region selector named %s in c3po.json!\n",name);
             error().throw_error_one(FLERR,"\nFATAL ERROR\n");
            }
        QJsonObject  selectorData = input().mainDoc()->object()[name].toObject();
        
   if(selectorData["fieldForSelection"].isNull())
     error().throw_error_one(FLERR,"\nYou have to specify the field for selection in your region selector!\n");
   
   QString qsV=selectorData["fieldForSelection"].toString();
   nameOfFieldForSelection_ = qsV.toUtf8().constData();
   
   if(selectorData["threshold"].isNull())
    error().throw_error_one(FLERR,"\nYou have to specify the threshold in your region selector!\n");
   
   threshold_ = selectorData["threshold"].toDouble();
   
   if(!selectorData["identifyBelowThreshold"].isNull())
    identifyBelowThreshold_=selectorData["identifyBelowThreshold"].toBool();
    
   if(!selectorData["updateRegionsEachStep"].isNull())
    updateRegionsEachStep_=selectorData["updateRegionsEachStep"].toBool();
    
    if(identifyBelowThreshold_) check=&SelectorCellRegion::checkIfBelow;
    else check=&SelectorCellRegion::checkIfAbove;
}

SelectorCellRegion::~SelectorCellRegion()
{
}
/*---------------------------------------------------------------------------*/
void SelectorCellRegion::begin_of_step()
{
  fieldForSelection_=dataStorage().fSFid(nameOfFieldForSelection_);
  if(fieldForSelection_==-1)
   error().throw_error_one(FLERR,"\nSelector Region: Your selection field has to be registered in CPPPO!\n");
   
  updateRegions=&SelectorCellRegion::IJKRegionSelector;
   

}
/*--------------------------------------------------------------------------*/
void SelectorCellRegion::middle_of_step()
{
  if(updated_) return;
 
 NofBubbles_=0;
 
 bubble_index_store.clear();
 bubble_index_eval=0;
 bubble_index_surface.clear();
 
 AA = new int[mesh().NofCells()];
   for(int selCell_=0;selCell_<mesh().NofCells();selCell_++)
     AA[selCell_]=0;
 
 for(int selCell_=0;selCell_<mesh().NofCells();selCell_++)
 {
  
  if(AA[selCell_]==1) continue;
  
  
  if( !(this->*check)(selCell_)) continue;
 
  bubble_index_store.push_back(selCell_); 
   
 
  while(int(bubble_index_store.size()) > bubble_index_eval)
   (this->*updateRegions)();
  
  createRegion();
 }
 
 if(!updateRegionsEachStep_) updated_ = true;
 delete AA;
}
/*--------------------------------------------------------------------------*/
void SelectorCellRegion::IJKRegionSelector()
{
 
 int cellId_= bubble_index_store[bubble_index_eval];
 bubble_index_eval++;
 AA[cellId_] = 1;
 int IJKcell[3];
 int x,y,z,id; 
 

   //check neighbors
   CellID2ijk(cellId_, &IJKcell[0]);
   bool boundary=false;
   for(int i=-1;i<2;i++)
   {
    x = IJKcell[0] + i;
    for(int j=-1;j<2;j++)
    {
     y = IJKcell[1] + j;
     for(int k=-1;k<2;k++)
     { 
      z = IJKcell[2] + k;
      
      id = ijk2CellIDOF(x,y,z);
      if(id == -1)
      {
       boundary=true;
       continue;
      }
      
      if(AA[id]!=0) continue;
      
      if(!(this->*check)(id))
      {
       boundary=true;
       continue;
      }
      
      //checkIndexStore(id);
      bubble_index_store.push_back(id);
      AA[id]=2;    
     }
    }
   }
 
  if(boundary) bubble_index_surface.push_back(cellId_);
   
}
/*-------------------------------------------------------------------------*/
bool SelectorCellRegion::checkIfBelow(int cell)
{
 
  if(dataStorage().SF(fieldForSelection_,cell) < threshold_ + TOLERANCE)
   return true;
   
  return false;

}
/*-------------------------------------------------------------------------*/
bool SelectorCellRegion::checkIfAbove(int cell)
{
 
  if(dataStorage().SF(fieldForSelection_,cell) > threshold_ - TOLERANCE)
   return true;

  return false;

}

/*-------------------------------------------------------------------------*/
int SelectorCellRegion::ijk2CellIDOF(int x, int y, int z) 
{
    
    int * bound = input().getBound();
    int ijk[3]={x,y,z};
    
    for(int i=0;i<3;i++)
    {
     
     if( bound[i]==-1)
     {
      if( ijk[i] >=  mesh().cells_subdomain()[i]) return -1;
      if( ijk[i] < 0 ) return -1;
     }
     
     if( bound[i]==-2)
     {
      if( ijk[i] >=  mesh().cells_subdomain()[i]) ijk[i] -= mesh().cells_subdomain()[i];
      if( ijk[i] < 0 )  ijk[i] += mesh().cells_subdomain()[i];
     }
           
    }
    
    int CellID_OF;          
   
    CellID_OF =   (ijk[2])*(  mesh().cells_subdomain()[0]  * mesh().cells_subdomain()[1]   )
                + (ijk[1])* mesh().cells_subdomain()[0] 
                +  ijk[0];
    
    return CellID_OF;
}
/*-------------------------------------------------------------------------*/
void SelectorCellRegion::CellID2ijk(int center_Cell_ID, int* result)
{
 ////////////////////////   Cell ID to i-j-k indexing   ///////////////////////////
   
 int x=mesh().cells_subdomain()[0];  
    
 int y=mesh().cells_subdomain()[1];
 result[2] = int(center_Cell_ID / (x*y)); //estimation of z-layer
 result[1] = int((center_Cell_ID-(x*y*result[2]))/x); //estimation of y-layer
 result[0] = int((center_Cell_ID-(x*y*result[2]))-((result[1]*x))); //estimation of x-layer
}
/*-------------------------------------------------------------------------*/
void SelectorCellRegion::checkIndexStore(int id)
{
 
 for(unsigned int i=0;i<bubble_index_store.size();i++)
  if(id==bubble_index_store[i])
   return;
 
 bubble_index_store.push_back(id);
  
}
/*-------------------------------------------------------------------------*/
void SelectorCellRegion::createRegion()
{
 
 char buf[140];
 
 sprintf(buf,"%s_%i",name(),NofBubbles_);
 
 dataStorage().createRegion(buf);
 
 dataStorage().selectRegion(dataStorage().NofRegions()-1);
 
 for( int i=0;i<bubble_index_eval;i++)
  dataStorage().addCellToRegion(bubble_index_store[i],false);
  
 for(unsigned int i=0;i<bubble_index_surface.size();i++)
  dataStorage().addCellToRegion(bubble_index_surface[i],true);
 
 bubble_index_surface.clear();
 bubble_index_store.clear();
 bubble_index_eval=0;
 
 NofBubbles_++;
 

}
