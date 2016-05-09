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
	Basic Class for 1D Binning
-----------------------------------------------------------------------------------*/

#ifdef OPERATION_CLASS

OperationStyle(binning, OperationBinning)

#else


#ifndef C3PO_OPERATION_BINNING_H
#define C3PO_OPERATION_BINNING_H

#include "operation_base.h"
#include "running_stat.h"
#include "qjson_includes.h"


#define BINCOUNT 6   //Default bin count

namespace C3PO_NS
{


class OperationBinning : public OperationBase
{
    public:

      OperationBinning(c3po *ptr,const char *_name);
      ~OperationBinning();

      void init(int narg, char const* const* arg);
      void init(QJsonObject jsonObj) ;

      void begin_of_step() {processSamples();};
      void middle_of_step() {};
      void end_of_step();


      void insertSample( std::vector<double>* y, std::vector<double>* x);
      void processSamples();
      void selectRunningStat();
      
      bool isFree() {return !hasSample_;}; 

    private:
    
    class        runningStat     *rs_ ;  
    std::vector< runningStat * >  rsVec_ ; //class for running statistics, 

    mutable int*         bincount_;     //number of bins, must come from JSON input script
    mutable double*      binLowBorder_; //lower bin boarder, must come from JSON input script
    mutable double*      binUpBorder_;  //upper bin boarder, must come from JSON input script
    mutable double*      binDelta_;     //bin delta

    mutable bool         overwrite_; //
    
    std::vector<double>*   marker_1;   //current pointer to marker
    std::vector<double>*    sample;     //current pointer to sample 
    const char*             Name_;
    
    mutable std::string dir_;
    
    mutable bool hasSample_;
    mutable int NofMarkers_;
    
    int sumProdArray(int* array,int* array2, int index);
    int prodArray(int* array, int index);
    int NofRs_;
};
} //end c3po_NS

#endif
#endif
