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

#ifndef C3PO_OPERATION_BASE_H
#define C3PO_OPERATION_BASE_H

#include "c3po_base_accessible.h"
#include "string.h"
#include "qjson_includes.h"


namespace C3PO_NS
{

//numbering of operation type
enum{ FILTER,		//0
      SAMPLE,   	//1
      BIN	        //2
    };


class OperationBase : public c3poBaseAccessible
{
    public:

      OperationBase(c3po *ptr,const char *_name);
      ~OperationBase();

      virtual void init(int narg, char const* const* arg) {};
      virtual void init(QJsonObject jsonObj)  {};

      virtual void begin_of_step() {};
      virtual void middle_of_step() {};
      virtual void end_of_step() {};

      virtual void insertSample(double, double) {}; //general interface to push values to an operation
      virtual void insertSample(double, int) {}; //general interface to push values to an operation
      virtual void insertSample(std::string, std::vector<double>*,std::vector<double>*,std::vector<double>*) {};
      virtual void insertSample(std::string,std::vector<double>*) {};
      virtual void insertSample(std::vector<double>*,std::vector<double>*) {};
      
      virtual void flush() {};              //general interface to flush/clear/clean an operation
      virtual void setFile(std::string) {};
      virtual void set() {};
     

      virtual int    sampleCount() {return 1;};     
      virtual double sampleDelta() {return 1.0;}; 
      
      virtual void automaticBinCount(int newBinCount) {};
       
      virtual std::vector<std::string> returnVectorFTFNames() {std::vector<std::string> v_; return v_;};
      virtual std::vector<std::string> returnScalarFTFNames() {std::vector<std::string> v_; return v_;};       
       
      virtual std::vector<std::string> returnVectorFTFVarianceNames() {std::vector<std::string> v_; return v_;};
      virtual std::vector<std::string> returnScalarFTFVarianceNames() {std::vector<std::string> v_; return v_;};              
      //Access functions
      const char* name() {return name_;}
      const char* scope() {  return scope_; }
      virtual double value(){return 0.0;};
      int getID() {return operationID_;};
      
      virtual bool particleBased() {return false;};
      
      virtual bool isFree() {return false;};

    protected:

      char *name_;
      char *scope_; // scope/name for JSON file

      mutable int   operationID_; //id of operation

      std::string   typeName_; //used to generate the dir name
      bool          typeNameDirGenerated_; //used to generate the dir name
};

} //end c3po_NS

#endif

