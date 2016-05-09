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

#ifndef C3PO_OPERATION_SAMPLING_H
#define C3PO_OPERATION_SAMPLING_H

#include "operation_base.h"
#include "output.h"
#include "qjson_includes.h"
#include "formula.h"


namespace C3PO_NS
{


class OperationSampling : public OperationBase
{
    public:

      OperationSampling(c3po *ptr,const char *_name);
      ~OperationSampling();

      virtual void init(QJsonObject jsonObj);
      virtual void process_input(QJsonObject jsonObj) {};
      
      void begin_of_step() {};
      void middle_of_step() ;
      void end_of_step() ;
      
      
      void insertSample(double _sample,double _x1,unsigned int id_=0, int marker_=0);
      void insertSample(double _sample,double* _x1,unsigned int id_=0);
      
      void setFile(std::string file)  {fileToWrite_.assign(file);};

      void flushToBin();
      void flushToFile();

      int    sampleCount() {return sampleCount_;};     
      double sampleDelta() {return sampleDelta_;}; 
      
      void clearSamples();
      void createSampleVectors(unsigned int, unsigned int, unsigned int);
      void sampleSetSkip(unsigned int id);
      
      void registerInputFields(std::string samplesVF, std::string samplesSF,  std::string markers) const;
      
     // void calculateFormula() {};
      
      
      
      

    protected:

      std::vector< std::vector<double> >   samples_;       //main storage for scalar-valued samples
      std::vector< std::vector<double>* >  samplesX1_;     //main storage for marker (x1...marker to identify bin)
      std::vector< std::string >           samplesNames_;  //names to indicate vector or scalar
      std::vector< bool >                  samplesSkip_;   //bool to indicate skipping of sampling (e.g., to suppress output)
      long double                          flushID_;       //id to mark samples

      mutable int                          sampleCount_;   //count of samples to be drawn (PER CELL or PER PARTICLE!)
      mutable double                       sampleDelta_;   //delta between samples

      mutable bool                         save2Bin_; 
      mutable bool                         save2Disk_;
      
      mutable std::string                  fieldName_;
      
      mutable std::vector<std::string>     VFtoSample_;
      mutable std::vector<std::string>     SFtoSample_;

      mutable std::vector<std::string>     markers_;
      
  
      std::string                          fileToWrite_;
      const char*                          Name_;
      mutable bool                         lagrangian_;
      
      int                                  binOp_;
      bool                                 overwrite_; 
      bool                                 execFormula_;
      
      Formula *                            formula_;
      mutable int                          NofMarkers_;
      
      std::string                          probesName_;

      bool                                 normalizeSample_;
      
      int                                  component_;
      
      bool                                 selective_;
      
      double                               minimum_[3];
      
      double                               maximum_[3];
            

};

} //end c3po_NS

#endif
