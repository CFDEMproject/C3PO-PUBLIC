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


#include "operation_binning.h"
#include "operation_container.h"
#include "selector_container.h"
#include "error.h"
#include <math.h>  
#include "comm.h"
#include "input.h"
#include "output.h"
#include <sstream>

using namespace C3PO_NS;

/* ----------------------------------------------------------------------
   ModelBase Constructor
------------------------------------------------------------------------- */

OperationBinning::OperationBinning(c3po *ptr,const char *name) 
: 
OperationBase(ptr,name),
//bincount_(1),
//binLowBorder_(0),
//binUpBorder_(1),
overwrite_(false),
marker_1(NULL),
sample(NULL),
Name_(name),
hasSample_(false),
NofMarkers_(1)
{
    operationID_ = operationContainer().binCount(); //set ID, this is the 'old' sampleCount!
    char buf[150];
    typeName_ = "c3po_binning";
    typeNameDirGenerated_ = false;
    sprintf(buf,"%s/%s_proc%i",typeName_.c_str(), Name_,comm().me());
    dir_.assign(buf);
    
}

// * * * * * * * * * * * * * * * * * * * * * * * * *
OperationBinning::~OperationBinning()
{
    for(unsigned int i=0;i<rsVec_.size();i++)
        delete rsVec_[i];

    rsVec_.clear();
    
    delete bincount_;
    delete binLowBorder_;
    delete binUpBorder_;
    delete binDelta_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * *
void OperationBinning::init(int narg, char const* const* arg)
{

   output().generateDir(typeName_, typeNameDirGenerated_);

}
// * * * * * * * * * * * * * * * * * * * * * * * * *
void OperationBinning::init(QJsonObject jsonObj) 
{
   if(!jsonObj["NumOfMarkers"].isNull())
    NofMarkers_=jsonObj["NumOfMarkers"].toInt();
     
    bincount_ = new int[NofMarkers_];     
    binLowBorder_ = new double[NofMarkers_]; 
    binUpBorder_ = new double[NofMarkers_];  
    binDelta_= new double[NofMarkers_];     

   if(jsonObj["NumOfMarkers"].isNull())
   {
    if( !jsonObj["bincount"].isNull() )
    {
      bincount_[0]= (jsonObj["bincount"].toInt());
    }
    else error().throw_error_one(FLERR,"You must specify a valid number of bins! \n");

    if( !jsonObj["binUpBorder"].isNull() ) binUpBorder_[0] =  jsonObj["binUpBorder"].toDouble();
    else binUpBorder_[0] = 1.0;
    
    if( !jsonObj["binLowBorder"].isNull() ) binLowBorder_[0]=  jsonObj["binLowBorder"].toDouble();
    else binUpBorder_[0] = 0.0;
    
    
    binDelta_[0] = (binUpBorder_[0]-binLowBorder_[0])/double(bincount_[0]);
    
    if(binDelta_==0) error().throw_error_one(FLERR,"ERROR: OperationBinning::init() You can not set \"binUpBorder\" equal to \"binLowBorder\"!");
    
   }
   else
   {
   
    
    for(int i=0;i<NofMarkers_;i++)
    {
     std::string MarkNum_("Marker");
     
     std::stringstream ss_;
     ss_ << i;
     std::string str = ss_.str();
     
     MarkNum_.append(str);
     
    // std::cout << "\n" << MarkNum_;
     
     //if( !jsonObj[MarkNum_.c_str()].isNull())  error().throw_error_one(FLERR,"Cannot find Marker\n");
     
     QJsonObject json2=jsonObj[MarkNum_.c_str()].toObject();
     
     if( !json2["bincount"].isNull() )
     {
      bincount_[i]= (json2["bincount"].toInt());
     }
     else error().throw_error_one(FLERR,"You must specify a valid number of bins! \n");

     if( !json2["binUpBorder"].isNull() ) binUpBorder_[i] =  json2["binUpBorder"].toDouble();
     else binUpBorder_[i] = 1.0;

     if( !json2["binLowBorder"].isNull() ) binLowBorder_[i]=  json2["binLowBorder"].toDouble();
     else binLowBorder_[i] = 0.0;
  
     binDelta_[i] = (binUpBorder_[i]-binLowBorder_[i])/double(bincount_[i]);
    
     if(binDelta_==0) error().throw_error_one(FLERR,"ERROR: OperationBinning::init() You can not set \"binUpBorder\" equal to \"binLowBorder\"!");
    }
   }
   
   if( !jsonObj["overwrite"].isNull() )
        overwrite_ =  jsonObj["overwrite"].toBool();
  
}

// * * * * * * * * * * * * * * * * * * * * * * * * *
void OperationBinning::end_of_step()
{

    if(sample==NULL)
    {
     output().write_screen_one("WARNING: a binning operation is not connected to a sample operation!");
     return;
    }
 for(int k=0;k<NofRs_;k++)
 {
  rs_[k].setTime(dataStorage().getTimeName());
  rs_[k].computeGlobals(overwrite_);
  rs_[k].dumpBinCenters(binLowBorder_[k],binDelta_[k],overwrite_);  
 }
}

// * * * * * * * * * * * * * * * * * * * * * * * * *
void OperationBinning::insertSample(std::vector<double>* y, std::vector<double>* x)
{
    hasSample_=true;
    marker_1 = y;
    sample   = x;
   
   
   for(int i=0;i<NofMarkers_;i++)
   {
   //  std::cout << "\nMarkerSize id: " << i << " = "<< marker_1[i].size() << " Samplesize: " << sample->size() ;
     if( !(marker_1[i].size() == sample->size()) )
            error().throw_error_one(FLERR,"ERROR: marker_1->size() NOT EQUAL sample->size()");
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * 

void OperationBinning::processSamples()
{

    selectRunningStat();

    int size=sample->size();
    int binId;
   
   
    
    for(int it=0; it < size; it++) 
    { 
       int rsId=0;
        int coord_[NofMarkers_-1];
       //Calculate the bin id, i.e the position with respect to the first marker
        double  markerVal = marker_1[0][it];
        if(markerVal <= binLowBorder_[0]) binId=0;
        else if(markerVal >= binUpBorder_[0]) binId = bincount_[0] -1;
        else binId = floor ( ( markerVal - binLowBorder_[0] ) / binDelta_[0] );
     
     //calculate the rs id, i.e. the position in the parameter space (except for the first marker) 
     for(int i=1;i<NofMarkers_;i++)   
     {   
       double  markerVal = marker_1[i][it];
       if(markerVal <= binLowBorder_[i]) coord_[i-1]=0;
       else if(markerVal >= binUpBorder_[i]) coord_[i-1] = bincount_[i] -1;
       else coord_[i-1] = floor ( ( markerVal - binLowBorder_[i] ) / binDelta_[i] );
             
     }
     
     for(int i=0;i<(NofMarkers_-1);i++)
      rsId += coord_[i]*prodArray(&bincount_[1],i);
     //std::cout << "\nbinID: " << binId << " rsID: " << rsId << " bincount[0]: "<< bincount_[0] << " bincount[1]: "<<bincount_[1] << "markerVal: " << markerVal;
     rs_[rsId].add((*sample)[it],binId);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * *
int OperationBinning::sumProdArray(int* array,int* array2, int index)
{
  int sum=0;
  for(int i=0;i<index;i++)
   sum += (array[i])*(array2[i]);
   
  return sum;
}
int OperationBinning::prodArray(int* array, int index)
{
  int prod=1;
  for(int i=0;i<index;i++)
   prod = prod*array[i];
   
  return prod;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * *
void OperationBinning::selectRunningStat()
{
    unsigned int id_ = selectorContainer().getFilterid();
    if(id_<rsVec_.size())
    {
        rs_=rsVec_[id_];
    }
    else
    {
        std::string file_(dir_);
        std::string filtName_(selectorContainer().getFilterName());
        file_.append("_"+filtName_);
        NofRs_=1;
        for(int i=1;i<NofMarkers_;i++)
         NofRs_= NofRs_*bincount_[i];
         
        runningStat* r_ = new runningStat[NofRs_];
        
        if(NofRs_==1) r_[0].initialize(bincount_[0],input().dumpFormat().c_str(),file_.c_str(),Name_);
        else
        {
         int coord_[NofMarkers_-1];
         int binprod[NofMarkers_-1];
         for(int i=0;i<(NofMarkers_-1);i++)
         {
           coord_[i]=0;
           binprod[i]=prodArray(&bincount_[1],i);
         }
        
         
         for(int i=0;i<NofRs_;i++)
         {
          for(int k=0;k<(NofMarkers_-1);k++)
           coord_[k]=0;
          
          std::string newFile_(file_ + "_id_");
          for(int j=(NofMarkers_-2);j>-1;j--)
          {
           coord_[j]= lround( ( i - sumProdArray(&coord_[0],&binprod[0],NofMarkers_-1) )/(binprod[j]) ); 
           
           std::stringstream ss_;
           ss_ << coord_[j];
           std::string str = ss_.str();
           newFile_.append(str);
          }
          // cout << "\n Rs name : " << newFile_ << "  proc: " <<comm().me() << " binprod: " << binprod[0];;
           r_[i].initialize(bincount_[0],input().dumpFormat().c_str(),newFile_.c_str(),Name_);
         }
        }
        
      rsVec_.push_back(r_);
      rs_=r_;
    }
}
