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

#include "selector_cellIJK.h"
#include "operation_sampling.h"
#include "operation_container.h"
#include "data_storage.h"
#include "error.h"
#include "selector_container.h"
#include <fstream>
#include <sys/stat.h>
#include <iomanip>
#include <cmath>
#include "input.h"
#include "output.h"
#include <sstream>

#define PI 3.14159265359
#define SMALL_NUMBER 1e-10

#ifdef H5_LIB
#include "h5_c3po.h"
using namespace H5_C3PO_NS;
#endif

using namespace C3PO_NS;

/* ----------------------------------------------------------------------
   ModelBase Constructor
------------------------------------------------------------------------- */

OperationSampling::OperationSampling(c3po *ptr,const char *name) 
: 
OperationBase(ptr,name),
flushID_(0),
save2Bin_(false),
save2Disk_(false),
Name_(name),
lagrangian_(false),
binOp_(-1),
execFormula_(false),
NofMarkers_(1),
normalizeSample_(false),
component_(-1),
selective_(false)
{

    overwrite_ = true;

    operationID_ = operationContainer().sampleCount(); //set ID of corresponding sampling operations

    //run=&OperationSampling::sample; //set standard run mode

    typeName_ = "c3po_sampling";
    typeNameDirGenerated_ = false;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
OperationSampling::~OperationSampling()
{
 for(unsigned int i=0;i<samples_.size();i++)
     delete samplesX1_[i];
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
void OperationSampling::middle_of_step()
{

 if(execFormula_)
 {

    //compute and normalize new samples  
     std::vector<double> newsample_;
     std::vector< std::vector<double> > newmarkers_;
     for(int marker_=0;marker_<NofMarkers_;marker_++)
     {
      std::vector<double> tmp_;
      newmarkers_.push_back(tmp_);
     }
     
     formula_->evaluate(&samples_,&samplesX1_,&newsample_,&newmarkers_);

     for(int i=0;i<int(newsample_.size());i++)
       for(int marker_=0;marker_<NofMarkers_;marker_++)
        insertSample( newsample_[i],newmarkers_[marker_][i],samples_.size()-1,marker_); 

 }
 //perform normalization if necessary
 /*if(normalizeSample_)
 {
    formula_->normalize(&samples_, &samplesX1_);    
    if(input().verbose())
        output().write_screen_one("\nCPPPO normalized all your samples! ");
 }
 */
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
void OperationSampling::end_of_step()
{
 flushToBin();
 flushToFile();
 clearSamples();
}

/* ----------------------------------------------------------------------
   Init with QJson data
------------------------------------------------------------------------- */
void OperationSampling::init(QJsonObject jsonObj) 
{
  lagrangian_=jsonObj["lagrangian"].toBool();
  
  if(lagrangian_)
  {
   if(jsonObj["probesName"].isNull())
   error().throw_error_one(FLERR,"\nYou should specify a \"probesName\" for every \"lagrangian\" sampling operation");
   
   QString qsV=jsonObj["probesName"].toString();
   probesName_ =qsV.toUtf8().constData();
   
  

      
  }

   sampleCount_ = jsonObj["sampleCount"].toInt();
   sampleDelta_ = jsonObj["sampleDelta"].toDouble();
   overwrite_   = jsonObj["overwrite"].toBool();  
   
   if( !jsonObj["Formula"].isNull() )
   {
    execFormula_=true;
    if( !jsonObj["NormalizeFormula"].isNull() )
    {
     QJsonObject normalObj_ =  jsonObj["NormalizeFormula"].toObject(); 
     formula_= new Formula(jsonObj["Formula"].toString().toUtf8().constData(),normalObj_["Normalization"].toString().toUtf8().constData());
     normalizeSample_  = true;
     formula_->setMultBasicObject( &(basicMultiphaseQty()) );
	 
      
	   double referenceDp = normalObj_["dp"].toDouble();
	   double referenceRhoP = normalObj_["RhoP"].toDouble();
	   double referenceg = normalObj_["g"].toDouble();
	   double referenceetaf = normalObj_["etaFluid"].toDouble();
	   double referencerhof = normalObj_["rhoFluid"].toDouble();
	   int    dragLaw    = normalObj_["dragLaw"].toInt();

     basicMultiphaseQty().setupSettling(referenceDp, 
                                        referenceRhoP, 
                                        referenceg, 
                                        referenceetaf, 
                                        referencerhof, 
                                        dragLaw
                                       );
    
    }
    else
    {
     formula_= new Formula(jsonObj["Formula"].toString().toUtf8().constData());
     
    }
    
    if(input().verbose())
    {
     output().write_screen_one("\nCPPPO registered your formula : ");
     output().write_screen_one( (formula_->getFormula()).c_str());
    }
   
    bool absoluteNumerator_   = false;
    bool absoluteDenominator_ = false;
    
    
    if(!jsonObj["absoluteNumerator"].isNull())
      absoluteNumerator_=jsonObj["absoluteNumerator"].toBool();
      
    if(!jsonObj["absoluteDenominator"].isNull())
     absoluteDenominator_=jsonObj["absoluteDenominator"].toBool();
   
    formula_->setAbsoluteValue(absoluteNumerator_,absoluteDenominator_);
   
   }
   
 
  
  if(jsonObj["VFieldsToSample"].isNull() && jsonObj["SFieldsToSample"].isNull())
   error().throw_error_one(FLERR,
                                "ERROR: You must specify valid 'VFieldsToSample' or 'SFieldsToSample' for general sampling \n");
  if(jsonObj["marker"].isNull())
   error().throw_error_one(FLERR,
                                "ERROR: You must specify a valid 'marker' for general sampling \n");
  if(sampleCount_ < 0)
    output().write_screen_one("\nWARNING: Will save all Euler cells for 'general' sampling. This might be expensive. \n");

  if(!jsonObj["component"].isNull())
   component_ = jsonObj["component"].toInt();

 
  std::string sampleVF(jsonObj["VFieldsToSample"].toString().toUtf8().constData());
  std::string sampleSF(jsonObj["SFieldsToSample"].toString().toUtf8().constData());
  std::string markers(jsonObj["marker"].toString().toUtf8().constData());
  
  registerInputFields(sampleVF,sampleSF,markers);
  
  NofMarkers_=markers_.size();
         
  if(VFtoSample_.size()>0 && (component_>2 || component_==-1 ) )
   error().throw_error_one(FLERR,
                                "ERROR:Invalid vector component! Valid vector components are:\n 0 (= x )\n 1 (= y )\n 2 (= z ) \n");

  
  process_input(jsonObj);
  //set formula-based sampling
  if(execFormula_) 
  {
        createSampleVectors(VFtoSample_.size(),SFtoSample_.size(),1);

        //set limits
        if(!jsonObj["formulaLimiterNumerator"].isNull())
        {
            QJsonArray limitsTemp=jsonObj["formulaLimiterNumerator"].toArray();
            formula_->setLimiterNumerator( (limitsTemp.at(0)).toDouble(),  
                                           (limitsTemp.at(1)).toDouble()
                                         );
        }
        if(!jsonObj["formulaLimiterDenominator"].isNull())
        {
            QJsonArray limitsTemp=jsonObj["formulaLimiterDenominator"].toArray();
            formula_->setLimiterDenominator( (limitsTemp.at(0)).toDouble(),  
                                             (limitsTemp.at(1)).toDouble()
                                           );
        }
        if(!jsonObj["formulaLimiterSymmetry"].isNull())
        {
            QJsonArray symmetryTemp=jsonObj["formulaLimiterSymmetry"].toArray();
            formula_->setLimiterSymmetry( (symmetryTemp.at(0)).toBool(), (symmetryTemp.at(1)).toBool() );
            if(   (symmetryTemp.at(0)).toBool() 
               && ( (formula_->getLimitNumerator())[0]<0 || (formula_->getLimitNumerator())[1]<0 ) 
              )
               error().throw_error_one(FLERR,
                                       "ERROR: you cannot set use a formulaLimiterSymmetry[0]=true and have a negative limit for formulaLimiterNumerator! \n");
            if(   (symmetryTemp.at(1)).toBool() 
               && ( (formula_->getLimitDenominator())[0]<0 || (formula_->getLimitDenominator())[1]<0 ) 
              )
               error().throw_error_one(FLERR,
                                       "ERROR: you cannot set use a formulaLimiterSymmetry[1]=true and have a negative limit for formulaLimiterDenominator! \n");
        }

  }
  else createSampleVectors(VFtoSample_.size(),SFtoSample_.size(),0);  

  //Set skips for output
  if(!jsonObj["saveOnlyFormula"].isNull())
    if(jsonObj["saveOnlyFormula"].toBool())
    {
        if(!execFormula_)
           error().throw_error_one(FLERR,
                                "ERROR: you cannot set 'saveOnlyFormula' and not provide a formula! \n");
        for(uint iSample=0; iSample<(samples_.size()-1); iSample++)
          sampleSetSkip(iSample);
    }

  
  
  if(!jsonObj["selective"].isNull())
  {
   if(jsonObj["selective"].toBool())
     {
      selective_=true;
      if(jsonObj["max"].isNull())
       error().throw_error_one(FLERR,"You must specify the max domain size when using selective sampling operations. \n"); 
      
      if(jsonObj["min"].isNull())
       error().throw_error_one(FLERR,"You must specify the min domain size when using selective sampling operations. \n"); 
      
      for(int i=0;i<3;i++)
      {
       maximum_[i]=jsonObj["max"].toArray()[i].toDouble();
       minimum_[i]=jsonObj["min"].toArray()[i].toDouble();
      }
      
     }    
  
  }

   
   if(execFormula_)
    formula_->interpretFormula(VFtoSample_.size(),SFtoSample_.size());
  
   //One Sample vector is the default setting
   createSampleVectors(1,0,0);

 
    save2Disk_ = jsonObj["save2Disk"].toBool();
    if(save2Disk_)
    {
        output().write_screen_one("\nWARNING: Will save samples to disk. This might consume significant disk space. \n");
        char buf[140];
        sprintf(buf,"%s/%s_proc%i_",typeName_.c_str(),Name_,comm().me());
        fileToWrite_.assign(buf);
        output().generateDir(typeName_, typeNameDirGenerated_);
    }
    if( !jsonObj["save2Bin"].isNull() )
    {
        save2Bin_ = jsonObj["save2Bin"].toBool();
        if(save2Bin_)
            output().write_screen_one("Will save samples to binning array. Be sure you have initialized a corresponding binning operation. \n");
    }
 
 if(!lagrangian_)
 {
   //Tell dataStorage to make the required fields Global when they are available 
   for(unsigned int i=0; i<VFtoSample_.size();i++)
   if((input().getVFnumber(VFtoSample_[i]))==-1) dataStorage().addFieldToConvert(VFtoSample_[i]);
   
   for(unsigned int i=0; i<SFtoSample_.size();i++)
   if((input().getSFnumber(SFtoSample_[i]))==-1) dataStorage().addFieldToConvert(SFtoSample_[i]);
       
   for(unsigned int i=0; i<markers_.size();i++)
   if((input().getSFnumber(markers_[i]))==-1) dataStorage().addFieldToConvert(markers_[i]);
 }  
   
}
/*------------------------------------------------------------------------*/
void OperationSampling::registerInputFields(std::string samplesVF, std::string samplesSF,  std::string markers) const
{
  
  //register vector fields to sample
  for (unsigned int it=0; it<samplesVF.size();it++)
  {
        
        if (samplesVF[it] != ' ' )
        {  
          
          std::string name;
         
          while (samplesVF[it] != ' ')
            {
               name.push_back(samplesVF[it]);
               it++;  
               if (it==samplesVF.size()) break;           
               
            }  
          VFtoSample_.push_back(name);
         
        }
  }
  
  //register scalar fields to sample
   for (unsigned int it=0; it<samplesSF.size();it++)
  {
        
        if (samplesSF[it] != ' ' )
        {  
          
          std::string name;
         
          while (samplesSF[it] != ' ')
            {
               name.push_back(samplesSF[it]);
               it++;  
               if (it==samplesSF.size()) break;           
              
            }  
          SFtoSample_.push_back(name);
         
        }
  }
  
  //register markers (scalar fields)
   for (unsigned int it=0; it<markers.size();it++)
  {
        
        if (markers[it] != ' ' )
        {  
          
          std::string name;
         
          while (markers[it] != ' ')
            {
               name.push_back(markers[it]);
               it++;  
               if (it==markers.size()) break;           
              
            }  
          markers_.push_back(name);
         
        }
  }
  
  if(execFormula_)
   formula_->setMaxNumberOfFields( markers_.size(), SFtoSample_.size(),VFtoSample_.size());
  
}


/* ----------------------------------------------------------------------
   Sample Handling
------------------------------------------------------------------------- */
void OperationSampling::insertSample(double _sample, double _x1,unsigned int id_, int marker_) 
{
    samples_[id_].push_back  (_sample);
    samplesX1_[id_][marker_].push_back(_x1);
  
}
void OperationSampling::insertSample(double _sample, double* _x1,unsigned int id_) 
{
    samples_[id_].push_back  (_sample);
    
    for(int marker_=0;marker_<NofMarkers_;marker_++)
     samplesX1_[id_][marker_].push_back(_x1[marker_]);
  
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /
void OperationSampling::createSampleVectors(unsigned int nOfVectors,unsigned int nOfScalars, unsigned int nFormula)
{  
   int vecCounter=0;
   while(samples_.size()<nOfVectors)
   {
     std::vector<double> s_;
     std::vector<double> * x_ = new std::vector<double>[NofMarkers_];
     
     samples_.push_back(s_);
     samplesX1_.push_back(x_);
    
     samplesSkip_.push_back(false);
     std::string strA="vec";
     std::ostringstream Str;
     Str << vecCounter;
     strA.append(Str.str());
     samplesNames_.push_back(strA);
     vecCounter++;
   }

   int scalarCounter=0;
   while(samples_.size()<nOfScalars+nOfVectors)
   {
     std::vector<double> s_;
     std::vector<double> * x_ = new std::vector<double>[NofMarkers_];
     
     samples_.push_back(s_);
     samplesX1_.push_back(x_);
    
     samplesSkip_.push_back(false);
     std::string strA="scalar";
     std::ostringstream Str;
     Str << scalarCounter;
     strA.append(Str.str());
     samplesNames_.push_back(strA);
     scalarCounter++;
   }

   int formulaCounter=0;
   while(samples_.size()<nOfScalars+nOfVectors+nFormula)
   {
     std::vector<double> s_;
     std::vector<double> * x_ = new std::vector<double>[NofMarkers_];
     
     samples_.push_back(s_);
     samplesX1_.push_back(x_);
    
     samplesSkip_.push_back(false);
     std::string strA="formula";
     std::ostringstream Str;
     Str << formulaCounter;
     strA.append(Str.str());
     samplesNames_.push_back(strA);
     formulaCounter++;
   }

} 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /
void OperationSampling::sampleSetSkip(unsigned int id)
{
    samplesSkip_[id] = true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * /
void OperationSampling::clearSamples()
{
  for(unsigned int i=0;i<samples_.size();i++)
   {
    samples_[i].clear();
    for(int j=0;j<NofMarkers_;j++)
     samplesX1_[i][j].clear();
    
    }

} 
/* * * * * * * * * * * * * * * * * * * * * * * *  ** * * * * * * * * * * * * * */
void OperationSampling::flushToBin()  
{
  //Add to bins if bin exists
  if(!save2Bin_)
    return;

   int bincount_=operationContainer().binCount();
   if(bincount_<1)
        error().throw_error_one(FLERR,"You do not have a binning operation! Cannot flush samples. \n");  

   if(binOp_==-1)
    for(int op=0;op<bincount_;op++) 
     if( operationContainer().bin(op)->isFree())
     {
      binOp_=op;
      break;
     }
   
   if(binOp_==-1)  error().throw_error_one(FLERR,"Cannot find a free binning operation! \n");  
   int ssize_=samples_.size();
   int validSample = -1;
   for(int i = 0;i<ssize_;i++ ) 
   {
    if(samplesSkip_[i])
        continue;

    validSample++;

    if((binOp_+validSample)< bincount_)
    { 
      //printf("...binning (id = %d) and name %s \n", 
     //   operationID_, operationContainer().bin(operationID_)->name());
      
    //Must push marker first!
     operationContainer().bin(binOp_+validSample)->insertSample(samplesX1_[i], &samples_[i]);
     operationContainer().bin(binOp_+validSample)->begin_of_step();
    
    }
    else
     printf("WARNING: cannot find binning operation for sampling operation with ID %d. Samples are lost!! \n", 
              (operationID_+i));
  
  }   

}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void OperationSampling::flushToFile()  //For general samples
{
 if(!save2Disk_)
     return;

  int ssize_=samples_.size();
  
  
  for( int sId_=0;sId_<ssize_;sId_++)
  {

   if(samplesSkip_[sId_])
        continue;

   std::string localfile_(fileToWrite_);
   
   if(!overwrite_) 
   {
    localfile_.append("time");
    localfile_.append(dataStorage().getTimeName()+"_");
   }
   localfile_.append(selectorContainer().getFilterName());
   std::string datasetName_(name_);  
   localfile_.append("_");
   localfile_.append(samplesNames_[sId_]);
   
  #ifdef H5_LIB
   if(input().dumpFormat().compare("hdf5")==0)
   {

    localfile_.append(".h5");
    
   // cout << "\n you should see 8 of these messages!\n";
    
    //MPI_Barrier(MPI_COMM_WORLD);
    createH5file(localfile_); 
   
    int totdata=samples_[sId_].size();

    for(int i=0;i<NofMarkers_;i++)
    {
     char bufM[8];
     sprintf(bufM,"%i",i); 
     std::string Mid(bufM);
     std::string mark(datasetName_+"_marker"+Mid); 
     OneArrayToH5(localfile_,mark.c_str(),&samplesX1_[sId_][i][0],totdata); 
    }
     std::string val(datasetName_+"_values"); 
     OneArrayToH5(localfile_,val.c_str(),&samples_[sId_][0],totdata);
    
    
   }
  #endif
  if(input().dumpFormat().compare("json")==0)
  {
    const int totdata=samples_[sId_].size();
    double * data2 = new double[totdata];
    double ** data1 = new double*[NofMarkers_];
    
    for(int i=0;i<NofMarkers_;i++)
     data1[i] = new double[totdata];
    
    localfile_.append(".json");
    
    for (int i=0;i<totdata;i++)
    {
      for(int j=0;j<NofMarkers_;j++)   
       data1[j][i]= samplesX1_[sId_][j][i];
     
       data2[i]= samples_[sId_][i];
    } 
     
    
    std::vector<double*> datavec;
    for(int j=0;j<NofMarkers_;j++) 
     datavec.push_back(data1[j]);
    datavec.push_back(data2);
    
    std::vector<std::string> namevec;
    for(int j=0;j<NofMarkers_;j++) 
    {
     char bufM[8];
     sprintf(bufM,"%i",j); 
     std::string Mid(bufM);
     
     namevec.push_back("markers"+Mid);
    
    }
    namevec.push_back("samples");

    Output::createQJsonArrays(localfile_,datasetName_, namevec,datavec,totdata,overwrite_); 
    
    delete data1;
    delete data2; 
      
   } 

  }

}
