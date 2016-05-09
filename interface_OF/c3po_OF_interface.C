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

#include "c3po_OF_interface.H"
#include "fvMesh.H"
#include "polyMesh.H"
#include "fvCFD.H"
#include "c3po.h"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include <cmath>
#include <algorithm> 

#include "interpolationCellPoint.H"
#include "interpolationCell.H"

using namespace C3PO_NS;

/* ----------------------------------------------------------------------
   c3poOFInterface Constructors
------------------------------------------------------------------------- */
c3poOFInterface::c3poOFInterface
(
// const Foam::keyType& key,
    const Foam::fvMesh& mesh,
    MPI_Comm communicator
)
:
// key_(key),
    mesh_(mesh),
    myC3po_(NULL),
    communicator_(communicator),
    debug_(false),
    me_OF_(-1),
    nprocs_OF_(-1)
{	

    MPI_Comm_rank(communicator,&me_OF_);
    MPI_Comm_size(communicator,&nprocs_OF_);

    Info<< endl;    //Info is only called on CPU 0
    Info<< "***************************************" << endl;
    Info<< "***        C3PO OF Interface        ***" << endl;
    Info<< "*** (c) DCS GmbH, Linz, Austria     ***" << endl;
    Info<< "*** (c) TU Graz, Graz, Austria      ***" << endl;
    Info<< "***************************************" << endl;
    Info<< endl; 
    Info<< "...launching c3po-core now..." << endl;
    Info<< endl; 
    MPI_Barrier(communicator);
    
    myC3po_ = new c3po(0,NULL,communicator_);
    mCheck_=  new meshCheck(mesh_,myC3po_);

}

/* ----------------------------------------------------------------------
   c3poOFInterface Destructors
------------------------------------------------------------------------- */
c3poOFInterface::~c3poOFInterface()
{
    deleteParticles();
    delete myC3po_;
    delete mCheck_;
}

/* ---------------------------------------------------------------------- */
void c3poOFInterface::run()
{
 //gather field data and send to c3po
 registerAllFields();
 //run c3po core
 runC3po();

 //clean everything
 MPI_Barrier(MPI_COMM_WORLD);
 deleteC3POfields();
 
}

/* ---------------------------------------------------------------------- */
void c3poOFInterface::runIB()
{
 //gather field data and send to c3po
 registerAllFields();
 
 if(particleID_.size()>0)
  processParticles();
 //run c3po core
 runC3po();

 //clean everything
 MPI_Barrier(MPI_COMM_WORLD);
 deleteC3POfields();

 
}

/* ---------------------------------------------------------------------- */
void c3poOFInterface::registerAllFields()
{
 

 int numVF = myC3po_->getVFnamesNumber();
 for(int i=0;i<numVF;i++)
  registerVectorField(&(const_cast<volVectorField&>(mesh_.lookupObject<volVectorField> (myC3po_->getVFnames(i).c_str()))),
                      myC3po_->getVFnames(i).c_str()
                      );
  
 int numSF = myC3po_->getSFnamesNumber();
 for(int i=0;i<numSF;i++)
  registerScalarField(&(const_cast<volScalarField&>(mesh_.lookupObject<volScalarField> (myC3po_->getSFnames(i).c_str()))),
                      myC3po_->getSFnames(i).c_str()
                      );
  
  registerC3POfields();

}


/* ---------------------------------------------------------------------- */
void c3poOFInterface::createGradients()
{
  std::vector<std::string> GradScalList_ = myC3po_->getGradientScalarList();
  std::vector<std::string> GradVecList_  = myC3po_->getGradientVectorList();
  std::vector<std::string> shearList_  = myC3po_->getShearRateList();
  
  //Create Gradients of Scalar Fields
  for(unsigned int scal = 0; scal < GradScalList_.size(); scal++)
  {
   
   std::string name("grad");
   
   
   for(unsigned int sf=0; sf<nameScalarFieldsC3PO_.size();sf++)
   {
    
    //look for the field name
    if(GradScalList_[scal].compare(nameScalarFieldsC3PO_[sf])==0)
    {
     name.append(GradScalList_[scal].c_str());
    // scalarFieldsC3PO_[sf]->correctBoundaryConditions();
     volVectorField * gradField = new volVectorField(name.c_str(), fvc::grad( *scalarFieldsC3PO_[sf]) );
     vectorFieldsC3PO_.push_back(gradField);
     nameVectorFieldsC3PO_.push_back(name);
     myC3po_->registerVF(name.c_str(),&((*gradField)[0].component(0)),&((*gradField)[0].component(1)),&((*gradField)[0].component(2)),3); 
     
    }
   
   }
  
  }

  //Create Gradients of Vector Fields
  for(unsigned int vec = 0; vec < GradVecList_.size(); vec++)
  {
   
   std::string name("grad");
   
   for(unsigned int vf=0; vf<nameVectorFieldsC3PO_.size();vf++)
   {
   
    if(GradVecList_[vec].compare(nameVectorFieldsC3PO_[vf])==0)
    {
     
     name.append(GradVecList_[vec].c_str());
     
     volVectorField * baseField =  vectorFieldsC3PO_[vf];
     
     std::string coordName(name + "X");
     //baseField->correctBoundaryConditions();
     volVectorField * vecFieldx = new volVectorField(coordName.c_str(), fvc::grad(baseField->component(0)) );
     vectorFieldsC3PO_.push_back(vecFieldx);
     nameVectorFieldsC3PO_.push_back(coordName);
     myC3po_->registerVF(coordName.c_str(),&((*vecFieldx)[0].component(0)),&((*vecFieldx)[0].component(1)),&((*vecFieldx)[0].component(2)),3); 
     
     coordName.assign(name + "Y");
     volVectorField * vecFieldy = new volVectorField(coordName.c_str(), fvc::grad(baseField->component(1)) );
     vectorFieldsC3PO_.push_back(vecFieldy);
     nameVectorFieldsC3PO_.push_back(coordName);
     myC3po_->registerVF(coordName.c_str(),&((*vecFieldy)[0].component(0)),&((*vecFieldy)[0].component(1)),&((*vecFieldy)[0].component(2)),3); 
     
     
     coordName.assign(name + "Z");
     volVectorField * vecFieldz = new volVectorField(coordName.c_str(), fvc::grad(baseField->component(2)) );
     vectorFieldsC3PO_.push_back(vecFieldz);
     nameVectorFieldsC3PO_.push_back(coordName);
     myC3po_->registerVF(coordName.c_str(),&((*vecFieldz)[0].component(0)),&((*vecFieldz)[0].component(1)),&((*vecFieldz)[0].component(2)),3); 
     
    
    }
   } 
  }

  //Create shear rates from vector fields
  for(unsigned int sr = 0; sr < shearList_.size(); sr++)
  {
   
   std::string name("shearRate");
   
   for(unsigned int vf=0; vf<nameVectorFieldsC3PO_.size();vf++)
   {
   
    if(shearList_[sr].compare(nameVectorFieldsC3PO_[sr])==0)
    {
     
     name.append(GradVecList_[sr].c_str());
     
     volVectorField * baseField =  vectorFieldsC3PO_[sr];
     
     vectorField grads_[3];
     
     grads_[0] =  fvc::grad(baseField->component(0));
     grads_[1]  =  fvc::grad(baseField->component(1));
     grads_[2]  =  fvc::grad(baseField->component(2));
     
     volScalarField * shearStress_ = new volScalarField(name.c_str(),baseField->component(0));
     
     
     
     shearStress_->internalField()=0;
     
    
     forAll(mesh_.cells(),cellI)
     {
      double temp_=0; 
      for(int i=0;i<3;i++)
       for(int j=0;j<3;j++)
        {
         
         temp_ += grads_[i][cellI].component(j)*grads_[i][cellI].component(j);
        }

        temp_= std::pow(temp_,0.5);
        
        (*shearStress_)[cellI]=temp_;
      }
     
     scalarFieldsC3PO_.push_back(shearStress_);
     nameScalarFieldsC3PO_.push_back(name);
    
      myC3po_->registerSF(name,&((*shearStress_)[0]));
    
     
    
    }
   } 
  }


}

/* ---------------------------------------------------------------------- */
void c3poOFInterface::registerVectorField( volVectorField * field, const char * name)
{
  vectorFieldsIn_.push_back(field);
  nameVectorFieldsIn_.push_back(name);
}

/* ---------------------------------------------------------------------- */
void c3poOFInterface::registerScalarField( volScalarField * field, const char * name)
{
  scalarFieldsIn_.push_back(field);
  nameScalarFieldsIn_.push_back(name);
}

/* ---------------------------------------------------------------------- */
void c3poOFInterface::registerC3POfields()
{    
    unsigned int numVF = myC3po_->getVFnamesNumber();
    for (unsigned int i=0;i<numVF;i++)
     for (unsigned int n=0; n<vectorFieldsIn_.size();n++)
       if (nameVectorFieldsIn_[n]==myC3po_->getVFnames(i))  
        { myC3po_->GlobalVF(myC3po_->getVFnames(i),&((*vectorFieldsIn_[n])[0].component(0)),&((*vectorFieldsIn_[n])[0].component(1)),&((*vectorFieldsIn_[n])[0].component(2)),3);
     }   
   if(numVF!=vectorFieldsIn_.size()) FatalErrorIn("\n ERROR: You should register all Vector Fields in C3PO. \n"); 
       
       unsigned int numSF = myC3po_->getSFnamesNumber();
    for (unsigned int i=0;i<numSF;i++)
      for (unsigned int n=0; n<scalarFieldsIn_.size();n++)
       if (nameScalarFieldsIn_[n]==myC3po_->getSFnames(i))  
         {
          myC3po_->GlobalSF(myC3po_->getSFnames(i),&((*scalarFieldsIn_[n])[0]));
         }
  if(numSF!=scalarFieldsIn_.size()) FatalErrorIn("\n ERROR: You should register all Scalar Fields in C3PO. \n"); 
   
   myC3po_->setTime(mesh_.time().timeName());
} 


/* ---------------------------------------------------------------------- */
void c3poOFInterface::deleteC3POfields()
{    
  myC3po_->resetGlobalFields();
  
  nameVectorFieldsIn_.clear();
  vectorFieldsIn_.clear();
  
  nameScalarFieldsIn_.clear();
  scalarFieldsIn_.clear();
  
  
}


/* ---------------------------------------------------------------------- */
void c3poOFInterface::writeFields() 
{
 if(!myC3po_->writeFields()) return;
 
 for (unsigned int i=0;i<vectorFieldsC3PO_.size();i++)
 {
 // vectorFieldsC3PO_[i]->correctBoundaryConditions();
  vectorFieldsC3PO_[i]->write();
 }
 for (unsigned int i=0;i<scalarFieldsC3PO_.size();i++)
 {
 // scalarFieldsC3PO_[i]->correctBoundaryConditions();
  scalarFieldsC3PO_[i]->write(); 
 }
}

/* ---------------------------------------------------------------------- */
void c3poOFInterface::deleteFields() 
{
 myC3po_->resetFields();
 for (unsigned int i=0;i<vectorFieldsC3PO_.size();i++)
   delete vectorFieldsC3PO_[i];

  for (unsigned int i=0;i<scalarFieldsC3PO_.size();i++)
   delete scalarFieldsC3PO_[i];
   
  vectorFieldsC3PO_.clear();
  scalarFieldsC3PO_.clear(); 
  
  nameScalarFieldsC3PO_.clear();
  nameVectorFieldsC3PO_.clear();
       

}  

/* ---------------------------------------------------------------------- */
void c3poOFInterface::createFilteredFields(int id)
{
 std::string filterName_(myC3po_->getFilterName(id));
 int OpFilterNum_ = myC3po_->getOpFiltNum();
  for(int f_=0;f_<OpFilterNum_;f_++)
  {
   std::vector<std::string> vectorFields_(myC3po_->vectorFTF(f_));
   std::vector<std::string> scalarFields_(myC3po_->scalarFTF(f_));
   std::vector<std::string> vectorVariance(myC3po_->vectorFTFVariance(f_));
   std::vector<std::string> scalarVariance(myC3po_->scalarFTFVariance(f_));
   std::string OpName_(myC3po_->getOpFilterName(f_));
   for (unsigned int n=0; n<vectorFieldsIn_.size();n++) 
   {
    for (unsigned int i=0; i<vectorFields_.size();i++)
     if(vectorFields_[i].compare(nameVectorFieldsIn_[n])==0)
     {   
         std::string fieldName_(nameVectorFieldsIn_[n]); 
         fieldName_.append("_");
         fieldName_.append(OpName_.c_str()); 
         fieldName_.append("_");
         fieldName_.append(filterName_.c_str());      
         volVectorField *v_= new volVectorField( fieldName_.c_str(),*vectorFieldsIn_[n]);
        
         vectorFieldsC3PO_.push_back(v_);
         nameVectorFieldsC3PO_.push_back(fieldName_);
         myC3po_->registerVF(nameVectorFieldsIn_[n],&((*v_)[0].component(0)),&((*v_)[0].component(1)),&((*v_)[0].component(2)),3); 
      }
    
    for (unsigned int i=0; i<vectorVariance.size();i++)
     if(vectorVariance[i].compare(nameVectorFieldsIn_[n])==0)
     {   
         char buf[26];
         sprintf(buf,"_var%i_",i);
         std::string fieldName_(nameVectorFieldsIn_[n]); 
         fieldName_.append("_");
         fieldName_.append(OpName_.c_str()); 
         fieldName_.append(buf);
         fieldName_.append(filterName_.c_str());      
         volVectorField *v_= new volVectorField( fieldName_.c_str(),*vectorFieldsIn_[n]);
        
         vectorFieldsC3PO_.push_back(v_);
         nameVectorFieldsC3PO_.push_back(fieldName_);
         
         myC3po_->registerVF(nameVectorFieldsIn_[n],&((*v_)[0].component(0)),&((*v_)[0].component(1)),&((*v_)[0].component(2)),3); 
      }  
   }    
   
   for (unsigned int n=0; n<scalarFieldsIn_.size();n++)
   {
    for (unsigned int i=0; i<scalarFields_.size();i++)
     if(scalarFields_[i].compare(nameScalarFieldsIn_[n])==0)
     {
       
       std::string fieldName_(nameScalarFieldsIn_[n]);
       fieldName_.append("_");
       fieldName_.append(OpName_.c_str()); 
       fieldName_.append("_");
       fieldName_.append(filterName_.c_str());        
       volScalarField *s_=new volScalarField( fieldName_.c_str(),*scalarFieldsIn_[n]);
        
       scalarFieldsC3PO_.push_back(s_);
       nameScalarFieldsC3PO_.push_back(fieldName_);
         
       myC3po_->registerSF(nameScalarFieldsIn_[n],&((*s_)[0]));

     }  
    
    for (unsigned int i=0; i<scalarVariance.size();i++)
     if(scalarVariance[i].compare(nameScalarFieldsIn_[n])==0)
     {
       char buf[26];
         sprintf(buf,"_var%i_",i);
         
       
       std::string fieldName_(nameScalarFieldsIn_[n]);
       fieldName_.append("_");
       fieldName_.append(OpName_.c_str()); 
       fieldName_.append(buf);
       fieldName_.append(filterName_.c_str());        
       volScalarField *s_=new volScalarField( fieldName_.c_str(),*scalarFieldsIn_[n]);
        
       scalarFieldsC3PO_.push_back(s_);
       nameScalarFieldsC3PO_.push_back(fieldName_);
      
         
       myC3po_->registerSF(nameScalarFieldsIn_[n],&((*s_)[0]));
     }  
    
    
    }
  }
}

/* ---------------------------------------------------------------------- */
void c3poOFInterface::runFilter(int id)
{   
    createFilteredFields(id);
    myC3po_->runFilters(id);
    
    //Create Gradient fields if necessary
    createGradients();
}


/* ---------------------------------------------------------------------- */
void c3poOFInterface::runSampling() 
{ 
   myC3po_->runSampling();
}



/* ---------------------------------------------------------------------- */
void c3poOFInterface::runBinning() 
{
    myC3po_->runBinning();
}


/* ---------------------------------------------------------------------- */
void c3poOFInterface::runC3po() 
{
 int NofFilters_=myC3po_->numberOfFilters();
 myC3po_->preRunOperations();

 for(int id=0;id<NofFilters_;id++)
  {
   runFilter(id);
   writeFields();
   MPI_Barrier(MPI_COMM_WORLD);
   runSampling();
   runBinning();
   MPI_Barrier(MPI_COMM_WORLD);
   deleteFields();
  }
  
 myC3po_->postRunOperations();
} 

/* ---------------------------------------------------------------------- */
void c3poOFInterface::processParticles()
{
 
 //Very expensive function, shall be used with caution
 
 //Every processor must hold all particle positions and forces
 std::vector<double> posx_,posy_,posz_,forcex_,forcey_,forcez_;
 
 int displ[nprocs_OF_];
 int numfrags[nprocs_OF_];
 
 int size_=particleID_.size();
 
 MPI_Allgather(&size_,1,MPI_INT,numfrags,1,MPI_INT,MPI_COMM_WORLD);  
        
 int sum = 0;
     
   for (int i = 0; i < nprocs_OF_; ++i) {
         
        displ[i] = sum;
        sum +=numfrags[i]; 
        }
 for(int i=0;i<size_;i++)
 {
  posx_.push_back(parPos_[i][0]);
  posy_.push_back(parPos_[i][1]);
  posz_.push_back(parPos_[i][2]);
 }
 
 std::vector<double> bufx_(sum,0);
 std::vector<double> bufy_(sum,0);
 std::vector<double> bufz_(sum,0);
 
 
 MPI_Barrier(MPI_COMM_WORLD);  
 MPI_Allgatherv(&posx_[0], numfrags[me_OF_], MPI_DOUBLE, &bufx_[0], numfrags, displ, MPI_DOUBLE,MPI_COMM_WORLD);
 MPI_Allgatherv(&posy_[0], numfrags[me_OF_], MPI_DOUBLE, &bufy_[0], numfrags, displ, MPI_DOUBLE,MPI_COMM_WORLD);
 MPI_Allgatherv(&posz_[0], numfrags[me_OF_], MPI_DOUBLE, &bufz_[0], numfrags, displ, MPI_DOUBLE,MPI_COMM_WORLD);

 

   //Evaluate global forces

 for(unsigned int f=0;f<parForceTmp_[0].size();f++) 
 {
  std::vector<double> forcex_,forcey_,forcez_;
  for(int i=0;i<size_;i++)
  {
   forcex_.push_back(parForceTmp_[i][f][0]);
   forcey_.push_back(parForceTmp_[i][f][1]);
   forcez_.push_back(parForceTmp_[i][f][2]);
  }
  
  std::vector<double> fbufx_(sum,0);
  std::vector<double> fbufy_(sum,0);
  std::vector<double> fbufz_(sum,0);
  

   
 // MPI_Barrier(MPI_COMM_WORLD);  
  MPI_Allgatherv(&forcex_[0], numfrags[me_OF_], MPI_DOUBLE, &fbufx_[0], numfrags, displ, MPI_DOUBLE,MPI_COMM_WORLD);
  MPI_Allgatherv(&forcey_[0], numfrags[me_OF_], MPI_DOUBLE, &fbufy_[0], numfrags, displ, MPI_DOUBLE,MPI_COMM_WORLD);
  MPI_Allgatherv(&forcez_[0], numfrags[me_OF_], MPI_DOUBLE, &fbufz_[0], numfrags, displ, MPI_DOUBLE,MPI_COMM_WORLD);

 
  for(int i=0;i<size_;i++)
  {
 
   for(int j=0;j<sum;j++)
    if(j>=(displ[me_OF_]) && j<displ[me_OF_+1]) continue;
    else if(parPos_[i][0]>( bufx_[j] - 1e-09) && parPos_[i][0]<( bufx_[j] + 1e-09) )
    {
     if(parPos_[i][1]<( bufy_[j] - 1e-09) && parPos_[i][1] > (bufy_[j] + 1e-09)) continue;
     else if(parPos_[i][2]<(bufz_[j] - 1e-09) && parPos_[i][2]>(bufz_[j] + 1e-09)) continue;
     else
     {
      parForceTmp_[i][f][0]+=fbufx_[j];
      parForceTmp_[i][f][1]+=fbufy_[j];
      parForceTmp_[i][f][2]+=fbufz_[j]; 
     }
    }
  }
 }
 
 //Evaluate global scalars
 if(parScalarTmp_.size()>0)
  for(unsigned int s=0;s<parScalarTmp_[0].size();s++)
  {
   
   double commScal[size_];
   for(int i=0;i<size_;i++)
    commScal[i]=parScalarTmp_[i][s];
   
   double sbuf[sum];
   
   MPI_Allgatherv(commScal, numfrags[me_OF_], MPI_DOUBLE, sbuf, numfrags, displ, MPI_DOUBLE,MPI_COMM_WORLD);
  
    for(int i=0;i<size_;i++)
  {
 
   for(int j=0;j<sum;j++)
    if(j>=(displ[me_OF_]) && j<displ[me_OF_+1]) continue;
    else if(parPos_[i][0]>( bufx_[j] - 1e-09) && parPos_[i][0]<( bufx_[j] + 1e-09) )
    {
     if(parPos_[i][1]<( bufy_[j] - 1e-09) && parPos_[i][1] > (bufy_[j] + 1e-09)) continue;
     else if(parPos_[i][2]<(bufz_[j] - 1e-09) && parPos_[i][2]>(bufz_[j] + 1e-09)) continue;
     else
     {
      parScalarTmp_[i][s]+=sbuf[j];
     }
    }
  }
  
  } 
  
  
 
  
 for(int i=0;i<size_;i++)
 {
 //If the particle lies exactly on the boundary register it if it is the upper local boundary or the global lower boundary  
   if (((parPos_[i][0] <= mCheck_->min_Domain()[0]) && (parPos_[i][0] != mCheck_->min_Domain_Global()[0])) || (parPos_[i][0]> mCheck_->max_Domain()[0])) continue;
   
   if (((parPos_[i][1] <= mCheck_->min_Domain()[1]) && (parPos_[i][1] != mCheck_->min_Domain_Global()[1])) || (parPos_[i][1]> mCheck_->max_Domain()[1])) continue;
  
   if (((parPos_[i][2] <= mCheck_->min_Domain()[2]) && (parPos_[i][2] != mCheck_->min_Domain_Global()[2])) || (parPos_[i][2]> mCheck_->max_Domain()[2])) continue;
   

 
  myC3po_->addParticle(particlesName_, pard_[i], parPos_[i], &parVel_[i][0], &parForceTmp_[i],&parScalarTmp_[i] ,parTorque_[i]);
  

 }
  
}

/* ---------------------------------------------------------------------- */
void c3poOFInterface::registerParticleIB( std::string           groupName,  //Name of the corresponding group
                                          int                         id_,  //Particle id
                                          double                        m,  //Particle radius
                                          double*                    pos_,  //Particle position
                                          double*                    vel_,  //Particle velocity
                                          std::vector< double* >   force_,  //Vector containing particle forces
                                          std::vector< double >* scalars_,  //Vector containing particle scalars (e.g. interphase heat)
                                          double*                  torque   //Particle torque
                                        )
{
  std::vector<int>::iterator it;
  
  //Beware the ghost particles!
  for(it = particleID_.begin(); it != particleID_.end(); it++)
   if(id_==*it)
    return;
  
  
  particlesName_=groupName;
  particleID_.push_back(id_);
  parVel_.push_back(vel_);
  pard_.push_back(m);
  parPos_.push_back(pos_);
  parForceTmp_.push_back(force_);
 
  if(scalars_!=NULL) 
   parScalarTmp_.push_back(*scalars_);
  
  parTorque_.push_back(torque);  
 
}

/* ---------------------------------------------------------------------- */
void c3poOFInterface::registerParticle( std::string           groupName,  //Name of the corresponding group
                                        int                         id_,  //Particle id
                                        double                        m,  //Particle radius
                                        double*                    pos_,  //Particle position
                                        double*                    vel_,  //Particle velocity
                                        std::vector< double*> *  force_,  //Vector containing particle forces
                                        std::vector< double >* scalars_,  //Vector containing particle scalars (e.g. interphase heat)
                                        double*                  torque   //Particle torque
                                      )
{
 std::vector<int>::iterator it;
 
 //Beware the ghost particles!
 for(it = particleID_.begin(); it != particleID_.end(); it++)
  if(id_==*it)
   return;
   
 double TOLERANCE = 1e-10;
 
 //If the particle lies exactly on the boundary register it if it is the upper local boundary or the global lower boundary
 for (int i=0;i<3;i++)
  if (((pos_[i] <= mCheck_->min_Domain()[i] - TOLERANCE) && (pos_[i] != mCheck_->min_Domain_Global()[i])) || (pos_[i] > mCheck_->max_Domain()[i]+TOLERANCE))
   return;

 
 myC3po_->addParticle(groupName, m, pos_, vel_, force_,scalars_, torque);
 
}

/* ---------------------------------------------------------------------- */

void c3poOFInterface::deleteParticles()
{
 std::vector<int*>::iterator it;
 

 particleID_.clear();
 parForceTmp_.clear();
 parScalarTmp_.clear();
 parVel_.clear();
 pard_.clear();
 parPos_.clear();
 parTorque_.clear();
 
 myC3po_->deleteParticles();
}


