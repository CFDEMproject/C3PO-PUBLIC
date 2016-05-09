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
#include "probeStorage.h"
#include "input.h"
#include "mesh.h"
#include "data_storage.h"
#include "output.h"
#include "stdio.h"
#include "timer.h"

#ifdef H5_LIB
 #include "h5_c3po.h"
 using namespace H5_C3PO_NS;
#endif



using namespace C3PO_NS;

probeStorage::probeStorage(std::string groupName, std::vector<std::string> filterNames, bool readParticlesFromJson, c3po *ptr)
:
c3poBaseAccessible(ptr),
groupName_(groupName),
filterNames_(filterNames),
MaxNofPar_(0),
readParticlesFromJson_(readParticlesFromJson)
{

  nParticlesProc_=new int[comm().nprocs()];

}


probeStorage::~probeStorage()
{
 deleteParticles();   
 delete nParticlesProc_;
}

//--------------------------------------------------------------------------------------

/* ----------------------------------------------------------------------
   Add particles / cells to memory
------------------------------------------------------------------------- */
void probeStorage::addParticle(double d, double* pos, double* vel, std::vector < double* >* force,std::vector<double>* scalars, double* torque)
{
 Particle* par_ = new Particle();
 particles_.push_back(par_);
 
 par_->setradius(d);
 par_->setpos(pos);
 par_->setvel(vel);
 par_->setforce(force);
 par_->setId(particles_.size()-1);
 
 if(torque != NULL)
  par_->settorque(torque);
  
 if(scalars != NULL)
  par_->setscalar(scalars);
}
/* ---------------------------------------------------------------------------*/
void probeStorage::readParticles() const
{
 
 timer().stamp();

 if(!readParticlesFromJson_) return;
 
 input().readParticles(&parPos_, groupName_);
 
 for(unsigned int par=0;par<parPos_.size()/3;par++)
 {
 
  //Check if within processor boundaries
  if( ! (
         ( parPos_[3*par] < mesh().dLocalMax()[0]  )  
          &&
          
         ( parPos_[3*par] > mesh().dLocalMin()[0]  ) 
        )
    ) continue;
 
  if( ! (
         ( parPos_[3*par+1] < mesh().dLocalMax()[1]  )  
          &&
          
         ( parPos_[3*par+1] > mesh().dLocalMin()[1]  ) 
        )
    ) continue;

  if( ! (
         ( parPos_[3*par+2] < mesh().dLocalMax()[2]  )  
          &&
          
         ( parPos_[3*par+2] > mesh().dLocalMin()[2]  ) 
        )
    ) continue;

 
  Particle* par_ = new Particle();
  particles_.push_back(par_);
  
  par_->setpos(&parPos_[3*par]);
  par_->setId(par);
 
 }
 
// removeGhostParticles(); This function will badly affect parallel scalability!!
 timer().stamp(TIME_INPUT);
}

/* ---------------------------------------------------------------------------*/
void probeStorage::removeGhostParticles() const
{
 
 int idvec[particles_.size()];
 
 for(unsigned int i=0; i<particles_.size();i++)
  idvec[i]=particles_[i]->getId(); 
  
 int counts=0;
 int displs[comm().nprocs()];
 int recvcounts[comm().nprocs()]; 
 int size=particles_.size();
 
 MPI_Allgather(&size,1,MPI_INT,&recvcounts[0],1,MPI_INT,MPI_COMM_WORLD);
 
 for(int n=0; n<comm().nprocs(); n++)
 {
  displs[n]=counts;
  counts+=recvcounts[n];
 
 }

 int buf[counts];
 
 MPI_Allgatherv(&idvec[0], particles_.size(), MPI_INT, &buf[0], &recvcounts[0], &displs[0], MPI_INT, MPI_COMM_WORLD);

 for(int i=0;i<size;i++)
 {
  for(int n=0;n<comm().nprocs();n++)
  {
   if(n==comm().me()) continue;
   
   for(int s=displs[n];s<recvcounts[n]+displs[n];s++)
   {
   
    if(particles_[i]->getId() == buf[s])
     if(n<comm().me())
     {
      delete particles_[i];
      particles_.erase(particles_.begin()+i);
     
     }
   
   }
  
  }
 
 }
 
 //Now fill  MaxNofPar_ 
 
 MaxNofPar_=0;
 
 for(int n=0; n<comm().nprocs(); n++)
  if(MaxNofPar_<recvcounts[n])
    MaxNofPar_=recvcounts[n]; 
}
 
/* ---------------------------------------------------------------------------*/
void probeStorage::deleteParticles()
{
 for (unsigned int i=0;i<particles_.size();i++)
  delete particles_[i];
 particles_.clear();
}
/* ---------------------------------------------------------------------------*/
void probeStorage::gatherParticleData() const
{
 int localNPar_ = particles_.size();
 
 MPI_Allgather(&localNPar_,1,MPI_INT,nParticlesProc_,1,MPI_INT,MPI_COMM_WORLD);
 
 //find the max
 for(int i=0;i<comm().nprocs();i++)
  if(MaxNofPar_<nParticlesProc_[i]) MaxNofPar_=nParticlesProc_[i];
 
 
}
//------------------------------------------------------------------------------
void probeStorage::writeParticles(std::string fileName_)
{

  timer().stamp();
  
  std::string file_(fileName_);
  
  file_.append("_");
  
  file_.append(groupName_);

  int npar=particles_.size();
  if(npar==0) return;
  
  
  #ifdef H5_LIB 
  if(input().dumpFormat().compare("hdf5")==0)
  { 
   file_.append(".h5");
   createH5file(file_); 
  
  std::vector< double > data; 
  for (int par=0; par<npar; par++)
   data.push_back(particles_[par]->getId());
 
   OneArrayToH5(file_, "Id", &data[0], npar);
   data.clear();
  
  for (int par=0; par<npar; par++)
   data.push_back(*particles_[par]->getradius());
 
   OneArrayToH5(file_, "radius", &data[0], npar);
   data.clear();
  for (int par=0; par<npar; par++)
   data.push_back(particles_[par]->getpos()[0]);
   
   OneArrayToH5(file_, "pos_x", &data[0], npar);
   data.clear();
  for (int par=0; par<npar; par++)
   data.push_back(particles_[par]->getpos()[1]);
   
   OneArrayToH5(file_, "pos_y", &data[0], npar);
   data.clear();
  for (int par=0; par<npar; par++)
   data.push_back(particles_[par]->getpos()[2]);
   
   OneArrayToH5(file_, "pos_z", &data[0], npar);
   data.clear();
  for (int par=0; par<npar; par++)
   data.push_back(particles_[par]->getvel()[0]);
  
   OneArrayToH5(file_, "vel_x", &data[0], npar);
   data.clear();
  
  for (int par=0; par<npar; par++)
   data.push_back(particles_[par]->getvel()[1]);
   
   OneArrayToH5(file_, "vel_y", &data[0], npar);
   data.clear();
  
  for (int par=0; par<npar; par++)
   data.push_back(particles_[par]->getvel()[2]);
   
   OneArrayToH5(file_, "vel_z", &data[0], npar);
   data.clear();
   
  for (int par=0; par<npar; par++)
   data.push_back(particles_[par]->gettorque()[0]);
  
   OneArrayToH5(file_, "torque_theta", &data[0], npar);
   data.clear();
  
  for (int par=0; par<npar; par++)
   data.push_back(particles_[par]->gettorque()[1]);
   
   OneArrayToH5(file_, "torque_phi", &data[0], npar);
   data.clear();
  

  int NofForces_=particles_[0]->getNofForces();
  for(int force=0;force<NofForces_;force++) 
  {
   char buf[40];
   sprintf(buf,"force_%i",force);
   std::string dsetName_(buf);
   for (int par=0; par<npar; par++)
    data.push_back((particles_[par]->getforce(force))[0]);
 
   dsetName_.append("_x");
   OneArrayToH5(file_, dsetName_.c_str() , &data[0], npar);
   data.clear();
   for (int par=0; par<npar; par++)
    data.push_back((particles_[par]->getforce(force))[1]);
  
   dsetName_.assign(buf);
   dsetName_.append("_y");
   OneArrayToH5(file_,dsetName_.c_str() , &data[0], npar);
   data.clear();
   dsetName_.assign(buf);
   dsetName_.append("_z");
   for (int par=0; par<npar; par++)
    data.push_back((particles_[par]->getforce(force))[2]);
  
   OneArrayToH5(file_,dsetName_.c_str() , &data[0], npar);
   data.clear();
  }  
 }
 #endif
 
 if(input().dumpFormat().compare("json")==0)
 {
  
   file_.append(".json");
  
  std::vector<double*>     dataV_;
  std::vector<std::string> dataN_;
  
  
   double dataId_[npar];
   double datar_[npar];
   double datax_[npar];
   double datay_[npar];
   double dataz_[npar];
   double datavx_[npar];
   double datavy_[npar];
   double datavz_[npar];
   double dataTp_[npar];
   double dataTt_[npar];
   int NofForces_=particles_[0]->getNofForces();
   double dataFx_[NofForces_][npar];
   double dataFy_[NofForces_][npar];
   double dataFz_[NofForces_][npar];
   for (int par=0; par<npar; par++)
   {
    dataId_[par]=particles_[par]->getId();
    datar_[par]=*particles_[par]->getradius();
    datax_[par]=particles_[par]->getpos()[0];
    datay_[par]=particles_[par]->getpos()[1];
    dataz_[par]=particles_[par]->getpos()[2];
    datavx_[par]=particles_[par]->getvel()[0];
    datavy_[par]=particles_[par]->getvel()[1];
    datavz_[par]=particles_[par]->getvel()[2];
    dataTp_[par]=particles_[par]->gettorque()[0];
    dataTt_[par]=particles_[par]->gettorque()[1];
    
    
    for(int force=0;force<NofForces_;force++) 
    {
     dataFx_[force][par]=particles_[par]->getforce(force)[0];
     dataFy_[force][par]=particles_[par]->getforce(force)[1];
     dataFz_[force][par]=particles_[par]->getforce(force)[2];
    }
 
   }
  dataV_.push_back(&dataId_[0]);
  dataN_.push_back("Id");
  
  dataV_.push_back(&datar_[0]);
  dataN_.push_back("radius");
  
  dataV_.push_back(&datax_[0]);
  dataN_.push_back("pos_x");
  
  dataV_.push_back(&datay_[0]);
  dataN_.push_back("pos_y");
  
  dataV_.push_back(&dataz_[0]);
  dataN_.push_back("pos_z");
  
  dataV_.push_back(&datavx_[0]);
  dataN_.push_back("vel_x");
  
  dataV_.push_back(&datavy_[0]);
  dataN_.push_back("vel_y");
  
  dataV_.push_back(&datavz_[0]);
  dataN_.push_back("vel_z");
  
  dataV_.push_back(&dataTp_[0]);
  dataN_.push_back("Torque_phi");
  
  dataV_.push_back(&dataTt_[0]);
  dataN_.push_back("Torque_theta");
  
  for(int force=0;force<NofForces_;force++) 
  {
   char buf[40];
   sprintf(buf,"force_%i_x",force);
   
   dataV_.push_back(&dataFx_[force][0]);
   dataN_.push_back(buf);
   
   
   sprintf(buf,"force_%i_y",force);
   
   dataV_.push_back(&dataFy_[force][0]);
   dataN_.push_back(buf);
  
   
   sprintf(buf,"force_%i_z",force);
   
   dataV_.push_back(&dataFz_[force][0]);
   dataN_.push_back(buf);
  
   
  }
  
  output().createQJsonArrays(file_,"particles",dataN_,dataV_, npar,true);
  

 }

 timer().stamp(TIME_OUTPUT);
}
//---------------------------------------------------------------//
void probeStorage::writeParticleFields(std::string fileName_)
{
  
  
  timer().stamp();

  std::string file_(fileName_);
  
  file_.append("_");
  
  file_.append(groupName_);

  int npar=particles_.size();
  if(npar==0) return;
  int fVF_size= dataStorage().fVFsize();
  int fSF_size= dataStorage().fSFsize();
  
 
  #ifdef H5_LIB  
 if(input().dumpFormat().compare("hdf5")==0)
 { 
  file_.append(".h5");
  createH5file(file_); 
  
  std::vector< double > data; 
 
  for(int i=0;i<fVF_size;i++)
   for(int j=0;j<3;j++)
   {
    for (int par=0; par<npar; par++)
     data.push_back(particles_[par]->filteredVector(i)[j]);
    
    char buf[40];
    sprintf(buf,"%s_%i",dataStorage().fVF(i)->name().c_str(),j);
    std::string dataset(buf);
    OneArrayToH5(file_, dataset.c_str(), &data[0], npar);
    data.clear();
   }

   for(int i=0;i<fSF_size;i++)
   {
    for (int par=0; par<npar; par++)
     data.push_back(*(particles_[par]->filteredScalar(i)));
    
    OneArrayToH5(file_, dataStorage().fSF(i)->name().c_str(), &data[0], npar);
    data.clear();
   }
   
    for (int par=0; par<npar; par++)
     data.push_back(particles_[par]->getFilterVolume());
    
    OneArrayToH5(file_, "filterVolume", &data[0], npar);
    data.clear();
   
   
  
 } 
 #endif 
 

if(input().dumpFormat().compare("json")==0)
 {
  //  error().throw_error_one(FLERR,"output of particle data in JSON format not yet supported");
  file_.append(".json");
  
  std::vector<double*> dataV_;
  std::vector<std::string> dataN_;
 
  for(int i=0;i<fVF_size;i++)
   for(int j=0;j<3;j++)
   {
    double * data = new double[npar];

    for (int par=0; par<npar; par++)
     data[par]=(particles_[par]->filteredVector(i)[j]);
    
    char buf[40];
    sprintf(buf,"%s_%i",dataStorage().fVF(i)->name().c_str(),j);
    std::string dataset(buf);
    dataN_.push_back(buf);
    dataV_.push_back(data);
    
   }

   for(int i=0;i<fSF_size;i++)
   {
    double * data = new double[npar];
    for (int par=0; par<npar; par++)
     data[par]=*(particles_[par]->filteredScalar(i));
    
    dataN_.push_back(dataStorage().fSF(i)->name());
    dataV_.push_back(data);
    
    }
   
    double * data = new double[npar];
    for (int par=0; par<npar; par++)
     data[par]=particles_[par]->getFilterVolume();
    
    dataN_.push_back("filterVolume");
    dataV_.push_back(data);
   
   output().createQJsonArrays(file_,"fields_at_particle_centers",dataN_,dataV_,npar,true);
   
   for(unsigned int i=0;i<dataV_.size();i++)
    delete dataV_[i];
    
   dataV_.clear();
 }
 
// cleanProbeFields();
 timer().stamp(TIME_OUTPUT);
}
//-----------------------------------------------------------------------------------------
bool probeStorage::runProbes(std::string filterName)
{
 if(filterNames_[0].compare("all")==0) return true;
 
 for(unsigned int i=0;i<filterNames_.size();i++)
  if(filterNames_[i].compare(filterName)==0) return true;
 
 return false;
}
//-------------------------------------------------------------------------------------------
void probeStorage::addScalar()
{
  for(unsigned int par=0;par<particles_.size();par++)
   particles_[par]->addScalar();

}

void probeStorage::addVector()
{
  for(unsigned int par=0;par<particles_.size();par++)
   particles_[par]->addVector();

}
//-------------------------------------------------------------------------------------------
void probeStorage::cleanProbeFields()
{
  for(unsigned int par=0;par<particles_.size();par++)
  {
   for(int vf=0;vf<particles_[par]->NofVF();vf++)
    for(int j=0;j<3;j++)
     particles_[par]->filteredVector(vf)[j]=0.0;
   
   for(int sf=0;sf<particles_[par]->NofSF();sf++)
     (*particles_[par]->filteredScalar(sf))=0.0;
 
     particles_[par]->setFilterVolume(0.0);
  
  }
}


