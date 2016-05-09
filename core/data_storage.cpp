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


#include "data_storage.h"
#include "input.h"
#include "comm.h"
#include "output.h"
#include "mesh.h"
#include "error.h"
#include "operation_container.h"
#include "selector_container.h"
#include <sstream>
#include <stdlib.h>
#include "timer.h"

#ifdef H5_LIB
 #include "h5_c3po.h"
 using namespace H5_C3PO_NS;
#endif


using namespace C3PO_NS;

/* ----------------------------------------------------------------------
  DataStorage Constructor
------------------------------------------------------------------------- */

DataStorage::DataStorage(c3po *ptr) :
   c3poBase(ptr),
      nproc(comm().nprocs()),
   timeName_("0"),
   currentProbe_(0),
   nbody_(0),
   nbody_all_(0),
   isAllocated_(false),
   useProbes_(false)
{
    haveFieldDir_ = false;
    haveParticleDir_ = false;
    haveRegionsDir_ = false;
    
    dirNameField_.assign(   "c3po_dataStorage_fields");
    dirNameParticle_.assign("c3po_dataStorage_particles");
    dirNameRegions_.assign( "c3po_dataStorage_regions");
   
} 

/* ---------------------------------------------------------------------------*/ 
DataStorage::~DataStorage()
{
    deleteRMA();
    deleteFields();
    deleteRegions();

    for(unsigned int probes=0;probes<probes_.size();probes++)
     delete probes_[probes];
     
    probes_.clear();
}

/////////////////////////////////////////////////////////////////////////////
                          // MEMBER functions
/////////////////////////////////////////////////////////////////////////////

/* ----------------------------------------------------------------------
   Intitialization phase - done before tun
------------------------------------------------------------------------- */

void DataStorage::init() const
{


 if(input().mainSettings()["useProbes"].isNull()) return;
 
 useProbes_ = input().mainSettings()["useProbes"].toBool();
 
 if(!useProbes_) return;
 
 QJsonDocument loadDoc;
 
 QJsonObject    parObj_;
 
 loadDoc = input().openJsonFile("c3po_control","probeSettings","probeSettings", parObj_ );
 
 //Read probe names
 
  std::vector<std::string> probeNames_;
 
  QString qsV=parObj_["probeNames"].toString();
  std::string big=qsV.toUtf8().constData();
   
  for (unsigned int it=0; it<big.size();it++)
  {
        
        if (big[it] != ' ' )
        {  
          
          std::string name;
         
          while (big[it] != ' ')
            {
               name.push_back(big[it]);
               it++;  
               if (it==big.size()) break;           
              
            }  
           
          probeNames_.push_back(name); 
           
        }
  }
  
  //Create probe storages
  
   for (unsigned int probe=0; probe<probeNames_.size();probe++)
   {
    
     if(parObj_[probeNames_[probe].c_str()].isNull())
      error().throw_error_one(FLERR,"\nOne or more defined probes are not specified in probeSettings.json\n");
      
     QJsonObject qOb=parObj_[probeNames_[probe].c_str()].toObject();
     
     //Read filters
     
     std::vector<std::string> filtNames_;
     
      if(qOb["filters"].isNull())
      error().throw_error_one(FLERR,"\nA \"filters\" entry is not specified in probeSettings.json\n");
 
     QString qsV1=qOb["filters"].toString();
     std::string big1=qsV1.toUtf8().constData();
   
     for (unsigned int it=0; it<big1.size();it++)
     {
        
        if (big1[it] != ' ' )
        {  
          
          std::string name;
         
          while (big1[it] != ' ')
            {
               name.push_back(big1[it]);
               it++;  
               if (it==big1.size()) break;           
              
            }  
           
          filtNames_.push_back(name); 
           
        }
     }
     
     
      if(qOb["readFromJson"].isNull())
      error().throw_error_one(FLERR,"\nA \"readFromJson\" entry is not specified in probeSettings.json\n");
 
    
     
     bool readFromJson = qOb["readFromJson"].toBool();
     
     //Create probe storage
     probeStorage * newStorage = new probeStorage( probeNames_[probe] ,
                                                   filtNames_,
                                                   readFromJson,
                                                   c3po_ptr()
                                                  );
                                                  
     probes_.push_back(newStorage);
        
   }

}


/* ----------------------------------------------------------------------
   Add particles / fields to memory
------------------------------------------------------------------------- */
void DataStorage::addParticle(std::string groupName,double d, double* pos, double* vel, std::vector < double* >* force, std::vector<double>* scalars , double* torque)
{
 if(!useProbes_) return;
 
 for(unsigned int probes=0;probes<probes_.size();probes++)
 {
  if(probes_[probes]->name().compare(groupName)==0)
  {
   probes_[probes]->addParticle(d,pos,vel, force, scalars,torque);
   return;
  }
 }
 
  error().throw_error_one(FLERR,"\nERROR: Invalid group name for particles/probes registration.\n");
 
}

/* ---------------------------------------------------------------------------*/
void DataStorage::deleteParticles()
{
 for (unsigned int i=0;i<probes_.size();i++)
  probes_[i]->deleteParticles();
}

/* ---------------------------------------------------------------------------*/ 
void DataStorage::writeFields(std::string OpName)
{

  if(!input().storageWriteFields())
      return;

  setFileName(OpName);
  output().generateDir(dirNameField_, haveFieldDir_);

  timer().stamp();
   
  int NX=mesh().NofCells();
  int fVF_size= fVF_.size();
  int fSF_size= fSF_.size();
  
 #ifdef H5_LIB 
 if( input().dumpFormat().compare("hdf5")==0 )
 {

  createH5file(filename_);
 
  for (int n=0;n<fVF_size;n++)    
  {
   double data[NX][3];
    
    for (int j = 0; j < NX; j++) 
	 for (int i = 0; i < 3; i++) 
	  data[j][i]=*(fVF_[n]->value(i,j));       
   std::string name(fVF_[n]->name());
   
   ThreeArrayToH5(filename_, name.c_str() , data, NX);
   
  }
  
  for (int n=0;n<fSF_size;n++)    
  {
   double data[NX];
    
    for (int j = 0; j < NX; j++) 
      data[j]=fSF_[n]->value()[j];       
   std::string name(fSF_[n]->name());
   
    OneArrayToH5(filename_, name.c_str() , data, NX);
   
  }

 }
 #endif
 if(input().dumpFormat().compare("json")==0)
 {

   std::vector<double**>    dataVPointer_;
   std::vector<std::string> dataNV_;
   
   std::vector<double*>     dataS_;
   std::vector<std::string> dataNS_;
   
   //save std::vector data
   for (int n=0;n<fVF_size;n++)     
   {
         dataVPointer_.push_back(fVF_[n]->values());
         std::string name(fVF_[n]->name());
         dataNV_.push_back(name);
   }
   std::string vname_(filename_);
   vname_.append("Vectors.json");
   if(dataNV_.size()>0)
     output().createQJsonArrays(vname_,dirNameField_,dataNV_, 
                              dataVPointer_, NX, 3, 
                              fVF_[0]->space()[0], true);

   //save scalar data
   for (int n=0;n<fSF_size;n++)    
   {
        std::string name(fSF_[n]->name());
        dataS_.push_back(fSF_[n]->value());
        dataNS_.push_back(name);   
   }
   std::string sname_(filename_);
   sname_.append("Scalars.json");
   if(dataNS_.size()>0)
     output().createQJsonArrays(sname_,dirNameField_,dataNS_,dataS_,NX,true);
    
 }

 timer().stamp(TIME_OUTPUT);
}


/* ---------------------------------------------------------------------------*/ 
void DataStorage::addfVF(std::string name,double* x,double* y,double* z, int sp)
{
 //create the field
  filteredVectorField *v_ = new filteredVectorField(name, x,y,z,sp);

 //reset values 
  for(int i=0; i < mesh().NofCells();i++)
   for(int j=0;j<3;j++)
    *v_->value(j,i)=0;
  
  fVF_.push_back(v_);
  
   for(unsigned int i=0;i<probes_.size();i++)
   {
    probes_[i]->addVector();
 
   }


}


/* ---------------------------------------------------------------------------*/ 
void DataStorage::addRMAvF(std::string name,double* x,double* y,double* z, int sp)
{
 
  filteredVectorField *v_ = new filteredVectorField(name, x,y,z,sp);

  RMAvF_.push_back(v_);
  
  
}

/* ---------------------------------------------------------------------------*/ 
void DataStorage::addfSF(std::string name,double* x)
{
 //create the field
  filteredScalarField *s_ = new filteredScalarField(name, x);
 //reset values
  for(int i=0; i < mesh().NofCells();i++)  
   s_->value()[i]=0;
  
  fSF_.push_back(s_);
  
 
  for(unsigned int i=0;i<probes_.size();i++)
 {
  probes_[i]->addScalar();
 
 }


}

/* ---------------------------------------------------------------------------*/ 
void DataStorage::addRMAsF(std::string name,double* x)
{

  filteredScalarField *s_ = new filteredScalarField(name, x);
 
  RMAsF_.push_back(s_);
 
  
}

/* ---------------------------------------------------------------------------*/ 
filteredVectorField* DataStorage::fVF(std::string name)
{
  int fVF_size= fVF_.size();
  
  for (int i=0;i<fVF_size;i++)
  {
  
    if (fVF_[i]->name()==name)
     return fVF_[i];
 
  }
  
 char buf[100];
 sprintf(buf,"Can not find a Vector field named: %s inside C3PO, have you registered it?", name.c_str());
 error().throw_error_all("data_storage.cpp",-1,buf);
 return NULL;

}

/* ---------------------------------------------------------------------------*/ 
filteredScalarField* DataStorage::fSF(std::string name)
{
  int fSF_size= fSF_.size();
  
  for (int i=0;i<fSF_size;i++)
  {
  
    if (fSF_[i]->name()==name)
     return fSF_[i];
    
    
  }
  
  char buf[100];
  sprintf(buf,"Can not find a Scalar field named: %s inside C3PO, have you registered it?", name.c_str());
  error().throw_error_all("data_storage.cpp",-1,buf);
  return NULL;

}


/* ---------------------------------------------------------------------------*/ 
void DataStorage::deleteFields()
{
 

 int fVF_size= fVF_.size();
 for (int i=0;i<fVF_size;i++)
   delete fVF_[i];
   
 fVF_.clear();
 
 int fSF_size= fSF_.size();
 for (int i=0;i<fSF_size;i++)
   delete fSF_[i];
   
 fSF_.clear(); 

}

/* ---------------------------------------------------------------------------*/ 
void DataStorage::deleteRMA()
{
 
  for(unsigned int i=0; i<RMAvF_.size();i++)
   delete RMAvF_[i];
  
  RMAvF_.clear();
  

  for(unsigned int i=0; i<RMAsF_.size();i++)
   delete RMAsF_[i];
  
  RMAsF_.clear();

}

/* ---------------------------------------------------------------------------*/ 
std::string DataStorage::NameChanges(std::string OpName)
{
   std::ostringstream temp[3];
   
   for (int i=0;i<3;i++)
    temp[i] << selectorContainer().filterWidth()[i];
   
   std::ostringstream temp2;
   temp2 << comm().me();
   
   std::string result("_"+timeName_+"_"+OpName.c_str()+"_processor"+temp2.str());
 
   return result;
}

/* ---------------------------------------------------------------------------*/ 
void DataStorage::setFileName(std::string OpName)
{


 if(input().dumpFormat().compare("hdf5")==0)
    filename_.assign(dirNameField_+"/results"+NameChanges(OpName)+".h5");
 
 if(input().dumpFormat().compare("json")==0)
    filename_.assign(dirNameField_+"/results"+NameChanges(OpName));

}

/* ---------------------------------------------------------------------------*/ 
void DataStorage::convertFilterToRMA(std::string name_)
{
 //delete if already existent 
 for(unsigned int i=0;i<RMAvF_.size();i++)
  if(name_.compare(RMAvF_[i]->name())==0 )
   {
    delete RMAvF_[i];
    RMAvF_.erase(RMAvF_.begin()+i);   
   }
   
 for(unsigned int i=0;i<RMAsF_.size();i++)
  if(name_.compare(RMAsF_[i]->name())==0 )
   {
     delete RMAsF_[i];
     RMAsF_.erase(RMAsF_.begin()+i);  
   }
 
 //Create RMA field
 for(unsigned int i=0;i<fVF_.size();i++)
  if(name_.compare(fVF_[i]->name())==0)
  {
   addRMAvF(name_,fVF_[i]->value(0,0),fVF_[i]->value(1,0),fVF_[i]->value(2,0), *fVF_[i]->space());
   return;
  }
 for(unsigned int i=0;i<fSF_.size();i++)
  if(name_.compare(fSF_[i]->name())==0)
  {
   addRMAsF(name_,fSF_[i]->value());
   return;
  }
  
}

/* ---------------------------------------------------------------------------*/ 
void DataStorage::refreshRMAfields()
{
 for(unsigned int i=0; i<fieldsToConvert_.size(); i++)
  convertFilterToRMA(fieldsToConvert_[i]);
}

/* ---------------------------------------------------------------------------*/ 
void DataStorage::addFieldToConvert(std::string field_)
{
 //Check is the field is already added to list
 for(unsigned int i=0;i<fieldsToConvert_.size(); i++)
  if(field_.compare(fieldsToConvert_[i])==0) return;
 
 fieldsToConvert_.push_back(field_);
}
/*-----------------------------------------------------------------------------*/
void DataStorage::writeParticles()
{
 if(!input().storageWriteParticles())
      return;

 setFileName("particles");
 output().generateDir(dirNameParticle_, haveParticleDir_);
 
 std::string file_(dirNameParticle_+"/time");
 file_.append(NameChanges("particles"));
 
 for(unsigned int i=0;i<probes_.size();i++)
 {
  probes_[i]->writeParticles(file_);
 
 }

}

/*------------------------------------------------------------------------------*/
void DataStorage::writeParticleFields(std::string filtName)
{
   if(!input().storageWriteParticles())
      return;

  setFileName(filtName);
  output().generateDir(dirNameParticle_, haveParticleDir_);

  std::string file_(dirNameParticle_+"/time");
  file_.append(NameChanges(filtName));
  
 
 for(unsigned int i=0;i<probes_.size();i++)
 {
  if(!(probes_[i]->runProbes(filtName))) continue;
  probes_[i]->writeParticleFields(file_);
 
 }

}

/*-------------------------------------------------------------------------------*/
void DataStorage::readParticles() const
{
 for(unsigned int i=0;i<probes_.size();i++)
 {
  probes_[i]->readParticles();
  probes_[i]->gatherParticleData();
 
 }

}
/*-------------------------------------------------------------------------------*/
void DataStorage::setProbes(std::string probeName) const
{
 for(unsigned int i=0;i<probes_.size();i++)
 {
  if(probeName.compare(probes_[i]->name())==0)
     currentProbe_=i;
 
 }

}
/*-------------------------------------------------------------------------------*/
void DataStorage::createRegion(std::string regionName)
{

 std::vector<std::string> filtName_;
 filtName_.push_back("all");
 
 region * reg_ = new region(regionName,filtName_,false,c3po_ptr());

 regions_.push_back(reg_);
}
/*-------------------------------------------------------------------------------*/
void DataStorage::addCellToRegion(int & IDInternal, bool isOnSurface)
{
 regionPtr_->addCell(IDInternal, isOnSurface);
}
/*-------------------------------------------------------------------------------*/
void DataStorage::deleteRegions()
{
 for(unsigned int i=0; i<regions_.size();i++)
  delete regions_[i];
  
 regions_.clear();

}
/*-------------------------------------------------------------------------------*/
void DataStorage::writeRegions()
{
 output().generateDir(dirNameRegions_, haveRegionsDir_);
 for(unsigned int i=0; i<regions_.size();i++)
  regions_[i]->write();
  
}
