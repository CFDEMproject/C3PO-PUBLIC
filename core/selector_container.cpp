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
#include "input.h"
#include "style_selector.h"
#include "error.h"
#include <vector>
#include <cmath>

using namespace C3PO_NS;

/* ----------------------------------------------------------------------
   SelectorContainer Constructor
------------------------------------------------------------------------- */

SelectorContainer::SelectorContainer(c3po *ptr) : c3poBase(ptr),
    selector_map_(new std::map<std::string,SelectorCreator>()),
    ijk(false),
    particleBased_(false)
{
  // fill map with all Selectors listed in style_selector.h

#define SELECTOR_CLASS
#define SelectorStyle(key,Class) \
  (*selector_map_)[#key] = &selector_creator<Class>;
#include "style_selector.h"
#undef SelectorStyle
#undef SELECTOR_CLASS


    filterWidth_[0]=-1;
    filterWidth_[1]=-1;
    filterWidth_[2]=-1;
    
 
    key_=new int[comm().nprocs()];

}

/* ---------------------------------------------------------------------- */

SelectorContainer::~SelectorContainer()
{
    delete selector_map_;
  
    delete key_;
   

}

/* ----------------------------------------------------------------------
   one instance per Selector in style_Selector.h
------------------------------------------------------------------------- */

template <typename T>
SelectorBase *SelectorContainer::selector_creator(c3po *ptr, char *name)
{
  return new T(ptr,name);
}

/* ----------------------------------------------------------------------
   settings
------------------------------------------------------------------------- */

void SelectorContainer::parse_command(int narg,char const* const* arg)
{
    int m = strlen(arg[0]) + 1;   
    char *SelectorType = new char[m];
    strcpy(SelectorType,arg[0]);

    int n = strlen(arg[1]) + 1;   
    char *SelectorName = new char[n];
    strcpy(SelectorName,arg[1]);
    

    if (selector_map_->find(arg[0]) != selector_map_->end())
    {
        SelectorCreator selector_creator = (*selector_map_)[arg[0]];
        if(strstr(SelectorType, "cell") != NULL)
        {
            cellSelector_.push_back(selector_creator(c3po_ptr(), SelectorName ) );
            cellSelector_[cellSelector_.size()-1]->init(narg, arg);
            
             if(strstr(SelectorType, "IJK") != NULL) ijk=true;
        }

        else if(strstr(SelectorType, "particle") != NULL)
        {
            particleSelector_.push_back(selector_creator(c3po_ptr(), SelectorName ) );
            particleSelector_[particleSelector_.size()-1]->init(narg, arg);
            
        }
        else if(strstr(SelectorType, "bubble") != NULL)
        {
            bubbleSelector_.push_back(selector_creator(c3po_ptr(), SelectorName ) );
            bubbleSelector_[bubbleSelector_.size()-1]->init(narg, arg);
            ijk=true;
        }
        else
        printf("FAIL: SelectorContainer PARSING: Selector type not properly set. Use 'cell', or 'particle' as the type. \n");
    }
    else if(strstr(SelectorType, "filter") != NULL)
    {
       //  printf("...registering filters sizes from file system/c3po.json\n");
         jsonDoc_ = input().mainDoc();
         if( jsonDoc_->isNull() )
             error().throw_error_one(FLERR,"QJsonDocument is invalid. Check! \n");
         
         if(jsonDoc_->object()[SelectorName].isNull())
            {
             printf("\nERROR: can not find the selector named %s in c3po.json!\n", SelectorName);
             error().throw_error_one(FLERR,"\nFATAL ERROR\n");
            }
         
         QJsonObject  filterData = jsonDoc_->object()[SelectorName].toObject();
    	 if( filterData.size() > 0 )
    	 {
    		filters_.push_back(new FilterBase(c3po_ptr(),SelectorName,-1) ); //Do not specify type here
            filters_[filters_.size()-1]->init(filterData);
    	 }
	     else
        	error().throw_error_one(FLERR,"Please add this entry in system/c3po.json:", SelectorName);
         
    }
    else
        printf("SelectorContainer PARSING: Selector name not found\n");
    
}


// ----------------------------------------------------------------------
void SelectorContainer::begin_of_step()
{
   cellSelector_[iOp]->begin_of_step();
  
}

// ----------------------------------------------------------------------
void SelectorContainer::middle_of_step()
{
   cellSelector_[iOp]->middle_of_step();

}
// ----------------------------------------------------------------------
void SelectorContainer::end_of_step()
{
  cellSelector_[iOp]->end_of_step();

}

// ----------------------------------------------------------------------
void SelectorContainer::bubble_run()
{
  
  for(unsigned int i=0;i<bubbleSelector_.size();i++)
  {
   bubbleSelector_[i]->begin_of_step();
   bubbleSelector_[i]->middle_of_step();
   bubbleSelector_[i]->end_of_step();
  }
}

// ----------------------------------------------------------------------
SelectorBase* SelectorContainer::selector(int i) const
{
    int cellSelector_size=cellSelector_.size();
    int cellSelector_count=selectorCount();
    
    
    if(i<cellSelector_size)
        return cellSelector_[i];

    else if(i < cellSelector_count )
        return particleSelector_[i-cellSelector_size];

    else
        return NULL;

}
// ------------------------------------------------------------------------
void SelectorContainer::checkIJK( bool struct_) const
{
 if(ijk==true && struct_==false)
  error().throw_error_one(FLERR, "\nERROR: your mesh is not structured but you choose a structured selector! Please use selector 'CellUnstruct' instead of 'cellIJK'.");

}
