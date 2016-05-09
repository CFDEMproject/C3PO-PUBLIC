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

#ifndef C3PO_INPUT_H
#define C3PO_INPUT_H


#include "stdio.h"
#include "error.h"
#include "input_base.h"
#include "c3po_base.h"
#include "c3po_base_interface.h"
#include "qjson_includes.h"
#include <fstream>
#include <string>


using std::ifstream;

namespace C3PO_NS
{

class Input : public c3poBase, public c3poBaseInterface, public InputBase
{
    public:

      Input(c3po *ptr);
      ~Input();

      void process_input_script();

      QJsonDocument openJsonFile(const char*, const char*, const char*, QJsonObject &json ) const;
      
      //return main QJsonDocument
      QJsonDocument *mainDoc()          const {return &loadDoc_;}
      QJsonObject   &mainSettings()     const {return mainObj_;}
      QJsonObject   &meshMainSettings() const {return meshMainSettings_;}
      
      bool verbose() {return mainObj_["verbose"].toBool();};
      bool meshverbose() {return meshMainSettings_["verbose"].toBool();};
      
      
      std::string getVFname(int x);
      std::string getSFname(int x);
      
      std::string dumpFormat() const {return dumpFormat_;};

      bool storageWriteFields()    const {return storageWriteFields_;};
      bool storageWriteParticles() const {return storageWriteParticles_;};
      
      
      int   getnumber_ofVF() {return VF_names.size();};
      int   getnumber_ofSF() {return SF_names.size();};
     
            
      void readFieldNames();
      
      int getVFnumber(std::string name) const
       { 
        for (unsigned int i=0; i<VF_names.size();i++)
          if (VF_names[i].compare(name)==0) return i;
          
         return -1;
       };
       
      int getSFnumber(std::string name) const
       { 
        for (unsigned int i=0; i<SF_names.size();i++)
          if (SF_names[i].compare(name)==0) return i;
        
        return -1;
       };
   
   
    double getMaxDomainfromJson(int i) const;
    double getMinDomainfromJson(int i) const;
    
    
    double cellSizefromJson(int i) const;
    
    inline int* getBound() const {return &bound_[0];};
    
    int getDomainDecomposition() const;
    
    bool writeInterface() const;
    
    void csv_columnsToRead( int *columns) const;   
    
    void readParticles(std::vector<double>* positions_, std::string probesName) const;
  
    bool registerCFDEMparticles() const;
    
    std::vector<std::string> getGradientScalarList() const {return SF_grad;};

    std::vector<std::string> getGradientVectorList() const {return VF_grad;};
    
    std::vector<std::string> getShearRateList() const {return shearRates_;};
 
    private:

     // istream infile_; 
      
     // std::ifstream ifnfile;                  // infile

      mutable QJsonDocument  loadDoc_;
      
      mutable QJsonObject    mainObj_;

      mutable QJsonObject    meshMainSettings_;
      
      

      mutable QString        myName;

      mutable int            writePrecision_;

      const char*            scope_;
      
      mutable int    bound_[3];
 
      mutable std::string    dumpFormat_;

      mutable bool storageWriteFields_;
      mutable bool storageWriteParticles_;
      
      mutable std::vector< std::string>  VF_names; //names of the std::vector fields
      mutable std::vector< std::string>  SF_names; //names of the scalar fields
      
      mutable std::vector< std::string>  VF_grad; //names of the gradient vector fields
      mutable std::vector< std::string>  SF_grad; //names of the gradient scalar fields
      
      mutable std::vector< std::string>  shearRates_;
      
      
      void readBCfromJson() const;
      
      
      

};

  // *************************************
  #include "input_I.h"
  // *************************************

} //end c3po_NS

#endif
