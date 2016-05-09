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
\*-----------------------------------------------------------------------------
    Description
    1-dimensional container for statistics
\*-----------------------------------------------------------------------------*/
#ifndef C3PO_RUNNING_STAT_H
#define C3PO_RUNNING_STAT_H

#include "stdio.h"
#include <string>
#include "mpi.h"

namespace C3PO_NS
{

class runningStat 
{
    public: 
        
        runningStat();
        ~runningStat();
        
        void    initialize(int size, const char* dumpFormat_, const char * directory, const char* OpName);
        
        void    add(double x, int bin);
        void    allocateMem(int _size);    
        void    updateFiles(bool);
        void    clear();

        double  *mean();
        double  *variance();

        int     *count();
        
        void     computeGlobals(bool overwrite);
        void     dumpBinCenters(double binLow, double delta, bool overwrite);
        
        void     setTime(std::string newTime_) const {time_.assign(newTime_);};
        
    private:
    
        int              size;          //number of means 
    
        //primary variables
        int            * count_;         //number of samples in each bin
        double         * run_mean_;      //mean in each bin
        double         * run_var_;       //variance in each bin
        double         * variance_;      //variance in each bin
    
        bool             willWriteToFile_;
        mutable int      writeIndex_;    //Index to identify writing step

        std::string file_;
        std::string file_g;

        mutable std::string dumpFormat;
        
        int me_, nprocs_;
        int num_;
        mutable std::string time_;
        mutable const char* OpName_;

        bool   binCentersDumped_;
        
        int      * bufCount;
        double   * bufMean;
        double   * bufVar;
};
}

#endif
