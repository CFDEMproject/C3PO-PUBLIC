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

#include "running_stat.h"
#include <cmath>
#include <sys/stat.h>
#include "memory_ns.h"
#include "output.h"
#include <string>

#ifdef H5_LIB
#include "h5_c3po.h"
using namespace H5_C3PO_NS;
#endif

using namespace C3PO_NS;
using namespace C3PO_MEMORY_NS;

runningStat::runningStat()
:
    size(1), //default size 
    count_(NULL),
    run_mean_(NULL),
    run_var_(NULL),
    num_(0)
{
   
   MPI_Comm_rank(MPI_COMM_WORLD,&me_);
   MPI_Comm_size(MPI_COMM_WORLD,&nprocs_);

   binCentersDumped_ = false;
   

} 

runningStat::~runningStat() 
{
    delete count_;
    delete run_mean_;
    delete run_var_;
    
    delete bufCount;
    delete bufMean;
    delete bufVar;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void runningStat::allocateMem(int _size)
{
    size = _size;
    create< int  >(count_,    size);
    create<double>(run_mean_, size);
    create<double>(run_var_,  size);
    create<double>(variance_, size);

    clear();
    
    create< int  >(bufCount,  nprocs_*size);
    create<double>(bufMean,   nprocs_*size);
    create<double>(bufVar,    nprocs_*size);
   


//    printf("running stats: mem allocated! \n");
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void runningStat::initialize(int size, const char* dumpFormat_,  const char * directory, const char* OpName)
{
     dumpFormat = dumpFormat_;
     OpName_ = OpName;
     file_.assign(directory);
     allocateMem(size);
   
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
void runningStat::computeGlobals(bool overwrite)
{
 
 double* var_=variance();

 std::string f(file_);
 
 //MPI_Barrier(MPI_COMM_WORLD);
 
 
 
 int displs[nprocs_];
 int recvCounts[nprocs_];
 
 for(int i=0;i<nprocs_;i++)
 {
  recvCounts[i]=size;
  displs[i]=i*size;
 }
   
   //Garher data in proc 0
   MPI_Gatherv(&count_[0], size, MPI_INT,
               &bufCount[0], recvCounts, displs,
                MPI_INT, 0, MPI_COMM_WORLD);
   
   MPI_Gatherv(&run_mean_[0], size, MPI_DOUBLE,
               &bufMean[0], recvCounts, displs,
                MPI_DOUBLE, 0, MPI_COMM_WORLD);
   
   MPI_Gatherv(&var_[0], size, MPI_DOUBLE,
               &bufVar[0], recvCounts, displs,
                MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
 
#ifdef H5_LIB
if(!dumpFormat.compare("hdf5"))
{
   
 double globalStat[size][3];
 
  //just for proc 0 
  if(me_==0)
  {
   //for every bin
   for(int i=0;i<size;i++)
   {
     globalStat[i][0]=0.;
     globalStat[i][1]=0.;
     globalStat[i][2]=0.;
    
    //sum processor values
    for(int p=0;p<nprocs_;p++)
    {
     
     globalStat[i][0]+=bufCount[i+p*size];                                           // Ctot= C1 + C2 + C3...
     globalStat[i][1]+=bufCount[i+p*size]*bufMean[i+p*size];                                    // Ctot*Mtot= C1*M1 + C2*M2 + C3*M3 +....
        
     globalStat[i][2]+=bufCount[i+p*size]*(bufVar[i+p*size] + bufMean[i+p*size]*bufMean[i+p*size]);     // Ctot*(Mtot^2 + Vtot) = C1(V1 + M1^2) + C2(... 
    }
    
    //check if != 0
    if(globalStat[i][0]>0)
     {
      globalStat[i][1]=globalStat[i][1]/globalStat[i][0];
      globalStat[i][2]=globalStat[i][2]/globalStat[i][0] - globalStat[i][1]*globalStat[i][1];
     }
   } 
  }

  
  
  if(me_==0) 
  {
     if(overwrite)    
     {
      f.append("_global.h5");
      createH5file(f);
     }
     else
     {
      char buf[40];
      sprintf(buf,"_global_time%s.h5",time_.c_str());
      f.append(buf);
      createH5file(f);
     }
     ThreeArrayToH5(f, OpName_, globalStat, size ); 
  
     }    
}  //end HDF5
 
#endif

if(!dumpFormat.compare("json"))
{
 
  std::vector<double*>      datavec;
  std::vector<std::string>  namevec;
 
  double globalStat[3][size];
 
  //just for proc 0 
  if(me_==0)
  {
   //for every bin
   for(int i=0;i<size;i++)
   {
     globalStat[0][i]=0.;
     globalStat[1][i]=0.;
     globalStat[2][i]=0.;
     
    //sum processor values
    for(int p=0;p<nprocs_;p++)
    {
     
     globalStat[0][i]+=bufCount[i+p*size];                                           // Ctot= C1 + C2 + C3...
     globalStat[1][i]+=bufCount[i+p*size]*bufMean[i+p*size];                                    // Ctot*Mtot= C1*M1 + C2*M2 + C3*M3 +....
        
     globalStat[2][i]+=bufCount[i+p*size]*(bufVar[i+p*size] + bufMean[i+p*size]*bufMean[i+p*size]);     // Ctot*(Mtot^2 + Vtot) = C1(V1 + M1^2) + C2(... 
    }
    
    //check if != 0
    if(globalStat[0][i]>0)
     {
      globalStat[1][i]=globalStat[1][i]/globalStat[0][i];
      globalStat[2][i]=globalStat[2][i]/globalStat[0][i] - globalStat[1][i]*globalStat[1][i];
     }
   } 
  }

  
  if(me_==0) 
  {
   if(overwrite)
    f.append("_global.json");
   else
   {
    char buf[40];
    sprintf(buf,"_global_time%s.json",time_.c_str());
    f.append(buf);
   } 
   
   datavec.push_back(globalStat[0]);
   namevec.push_back("count");
   datavec.push_back(globalStat[1]);
   namevec.push_back("mean");
   datavec.push_back(globalStat[2]);
   namevec.push_back("variance");
   
   Output::createQJsonArrays(f,OpName_, namevec,datavec, size ,overwrite); 
  }   
} //end JSON

 if(!overwrite)
  clear();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * *
void runningStat::dumpBinCenters(double binLow, double delta, bool overwrite)
{

    if( me_>0 || binCentersDumped_ )
        return;

    std::string f(file_);

    #ifdef H5_LIB
    if(!dumpFormat.compare("hdf5"))
    {
        double centers[size];
        
       
      if(overwrite)
       f.append("_binCenters.h5");
   
      else
      {
       char buf[40];
       sprintf(buf,"_binCenters_time%s.h5",time_.c_str());
       f.append(buf);
      } 
      
        for(int i=0;i<size;i++)
        { 
            centers[i] = binLow + (double(i)+0.5) * delta;
        }

        createH5file(f);
        OneArrayToH5(f, OpName_, centers, size ); 
    
    }  //end HDF5
    #endif

    if(!dumpFormat.compare("json"))
    {
        double centers[size];
        std::vector<double*>      datavec;
        std::vector<std::string>  namevec;
        
        if(overwrite)
         f.append("_binCenters.json");
   
        else
        {
         char buf[40];
         sprintf(buf,"_binCenters_time%s.json",time_.c_str());
         f.append(buf);
        } 
    
  
        for(int i=0;i<size;i++)
        { 
            centers[i] = binLow + (double(i)+0.5) * delta;
        }

        datavec.push_back(centers);
        namevec.push_back("centers");

        Output::createQJsonArrays(f,OpName_, namevec,datavec, size ,true); 
    }   //end JSON

    binCentersDumped_ = true;

} 

// * * * * * * * * * * * * * * * * * * * * * * * * * *
void runningStat::add(double x, int bin)
{
    if(bin>=size)
     return;
    count_[bin]    += 1;
    double old      = run_mean_[bin];
    run_mean_[bin] += (x - run_mean_[bin])/count_[bin];
    run_var_[bin]  += (x - old)*(x-run_mean_[bin]);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * *
void runningStat::updateFiles(bool overwrite)
{
 
  double* var_=variance();
  double data_[size][3];
   
  std::vector<double*>      datavec;
  std::vector<std::string>  namevec;
  
  std::string f(file_);
  
   #ifdef H5_LIB
    if(!dumpFormat.compare("hdf5"))
    {  
      if(overwrite)
        f.append(".h5");
   
      else
      {
       char buf[40];
       sprintf(buf,"_time%s.h5",time_.c_str());
       f.append(buf);
      } 
      
      createH5file(f);
    }
   #endif
   
   if(!dumpFormat.compare("json"))
   {
    
    if(overwrite)
      f.append(".json");
    else
    {
     char buf[40];
     sprintf(buf,"_time%s.json",time_.c_str());
     f.append(buf);
    } 

   for(int i=0;i<size;i++)
    { 
      data_[i][0]=count_[i];
      data_[i][1]=run_mean_[i];
      data_[i][2]=var_[i];
      
      if(!dumpFormat.compare("json"))
      {
       datavec.push_back(data_[i]);
       namevec.push_back("bin");
      }
    }

   #ifdef H5_LIB
   if(!dumpFormat.compare("hdf5"))
    ThreeArrayToH5(f, OpName_, data_, size );
   #endif
    if(!dumpFormat.compare("json"))
     Output::createQJsonArrays(f,OpName_, namevec,datavec,size,overwrite); 
     
   }
  
}

// * * * * * * * * * * * * * * * * * * * * * * * * * *
void runningStat::clear()
{
   for(int i=0;i<size;i++)
   {
        count_[i]    = 0;
        run_mean_[i] = 0.0;
        run_var_[i]  = 0.0;
   }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * *
int* runningStat::count() 
{

    return count_;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * *    
double* runningStat::mean()
{
    return run_mean_;
}    

// * * * * * * * * * * * * * * * * * * * * * * * * * *
double* runningStat::variance()
{

    for(int i=0;i<size;i++)
        variance_[i] =  (count_[i] > 1) ? run_var_[i]/(count_[i]-1) : 0.0 ;
    
    return variance_; 
}


