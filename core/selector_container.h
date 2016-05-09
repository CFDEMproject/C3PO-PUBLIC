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
	Container to hold selectors to determine the index of (i) Eulerian cells, and
    (ii) particle indices that we need to average over. Also holds a vector with the
    filters (i.e., the class that holds the size of the averaging region)
-----------------------------------------------------------------------------------*/

#ifndef C3PO_SELECTOR_CONTAINER_H
#define C3PO_SELECTOR_CONTAINER_H

#include "stdio.h"
#include "c3po_base.h"
#include "c3po_base_interface.h"
#include "filter_base.h"
#include "selector_base.h"
#include "output.h"
#include "comm.h"
#include <map>
#include <string>
#include "qjson_includes.h"

#include "data_storage.h"






namespace C3PO_NS
{

class SelectorContainer : public c3poBase, public c3poBaseInterface
{
    public:

      SelectorContainer(c3po *ptr);
      virtual ~SelectorContainer();
      
     

      void parse_command(int narg,char const* const* arg);

      void begin_of_step();
      void middle_of_step();
      void end_of_step();
      
      void bubble_run();
      
      int filterCount()   const {return     filters_.size();};
      int selectorCount() const {return     cellSelector_.size()
                                          + particleSelector_.size();
                                 };

      FilterBase*   filter(int idx)           {return filters_[idx];}
      SelectorBase* particleSelector(int idx) {return particleSelector_[idx];}
      SelectorBase* cellSelector(int idx)     {return cellSelector_[idx];}
      SelectorBase* bubbleSelector(int idx)     {return bubbleSelector_[idx];}

      SelectorBase* selector(int i)   const ;//return global Selector

      void setFilter(int x) const {
      
        iF=x;
        iOp=x;
        char message[256];
        sprintf(message,"\nRunning filter %s with ID: %d \n",
                 filters_[iF]->name(),
                 iF
               );
        output().write_screen_one(message);
      
        filters_[iF]->setNFilterCells();
        filterWidth_[0]=filters_[x]->filterNcells(0);
        filterWidth_[1]=filters_[x]->filterNcells(1);
        filterWidth_[2]=filters_[x]->filterNcells(2);

        
      
      };


      inline double* getFilterSize(int i) {return filters_[iF]->filterSize(i);};
      
      inline int* filterWidth() const {return &filterWidth_[0];};
      
      int *filtersize() const {  return currFilterCellCount_global_;};
      
      int filterType() const {return filters_[iF]->filterType();};
      
      bool selectiveFilter() const {return filters_[iF]->selectiveFilter(); };
      
      double* maxToFilter() const {return filters_[iF]->maxToFilter();};
      double* minToFilter() const {return filters_[iF]->minToFilter();};
   
      bool haveCell() {return cellSelector_[iOp]->haveCell();}	
      bool haveParticle() {return cellSelector_[iOp]->haveParticle();}	
      
      inline double*   currentCellCenters() {return cellSelector_[iOp]->currentCellCenters();}
             
	inline int* currentCell()           const { return &selectedCellID_;};
        int  currentLocCell()        const;
	void setCurrentCell(int i)          const {selectedCellID_=i;}; 
	
	int selectCell(double x, double y, double z)  const { return cellSelector_[iOp]->findNearestCell( x, y, z); };  
	
	void addCell(int i, int j)          const {(selected_Cells[j]).push_back(i);};
	void addFSize()                 const { };
	          
	inline int* getSCell(int i)         const {return &(selected_Cells[i][0]);};  
	inline double* getFSize(int i)         const {return &filterSize_[i];};
	
	void resetSCells()                  const {for (int i=0;i<comm().nprocs();i++) selected_Cells[i].clear();};
	void resetFSize()                  const {filterSize_.clear();};
	
	inline int SCsize(int i)            const {return selected_Cells[i].size();};
	
	std::string getFilterName() const {return filters_[iOp]->name();};
	
	int getFilterid() {return iOp;};
		
	void addCellInFilter(int id)          const {CellsInFilter_.push_back(id);};
	void writeKey(int proc, int value_)   const {key_[proc]=value_;};
	
	void resetCellsInFilter()               const {CellsInFilter_.clear();};
	
	inline double* filterVolume()  const {return &filterVolume_; };
	inline int filterVolumeSize()  const {return 1; };
	
	void addFilterVolume(double vol)  const { filterVolume_=vol; };
	void resetFilterVolume()          const {filterVolume_=0.;};
	
	inline std::vector<int>* getCellsInFilter() const {return &CellsInFilter_;};
	inline int*              getKey()           const {return key_;};    
	
	void                     setFSize(double size) const {filterSize_.push_back(size);};
	
	void checkIJK(bool struct_)   const;  
	
	void setCurrentParticle(int id_)  const {currentParticleId_=id_;} ;
	void isParticle()                 const {particleBased_=true;};
	void isCell()                     const {particleBased_=false;};
	
	int  currentParticle()          {return currentParticleId_;};   
	
	bool particleBased()            {return particleBased_;};
	                   
/*-----------------------------------------------------------------------------------------------------*/      
    void SCreport()
    {
       
    }	
/*-----------------------------------------------------------------------------------------------------*/
     private:

      typedef SelectorBase *(*SelectorCreator)(c3po*, char *name);
      std::map<std::string,SelectorCreator> *selector_map_;

      template <typename T> static SelectorBase *selector_creator(c3po *ptr, char *name);
      
     
      mutable double            filterVolume_;
           
      mutable std::vector<int> CellsInFilter_;
      mutable int*             key_;
      
     
      vector<FilterBase*>   filters_;
      vector<SelectorBase*> cellSelector_;
      vector<SelectorBase*> particleSelector_;
      vector<SelectorBase*> bubbleSelector_;
  
      mutable std::vector<int>* selected_Cells;
      mutable int      N_of_SCells;
      

      mutable std::vector<double> filterSize_;

      QJsonDocument *jsonDoc_;
      mutable int     selectedCellID_;
      mutable int filterWidth_[3];  //current filterWidth,  

      mutable int currFilterCellCount_global_[3];  //
 
      mutable int iF;
      mutable int iOp;
      
      mutable bool ijk;
      
      mutable int currentParticleId_;
      mutable bool particleBased_;

};

} //end c3po_NS

#endif
