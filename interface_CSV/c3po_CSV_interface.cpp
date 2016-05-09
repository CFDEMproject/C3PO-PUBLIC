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
#include "c3po_CSV_interface.h"
#include "qjson_includes.h"
#include <string>
#include <iostream>
#include <stdlib.h> 
#include <cmath>
#include <sstream>
#include <sys/stat.h>

#define MAX_CHARS_PER_LINE 1024

using namespace C3PO_NS;

struct stat dirLogicCSV;

/* ----------------------------------------------------------------------
   c3poCSVInterface Constructors
------------------------------------------------------------------------- */
c3poCSVinterface::c3poCSVinterface(MPI_Comm comm)
:
myC3po_(NULL),
comm_(comm),
comm_me_(0),
nprocs_(0),
csv_(NULL),
fields_(NULL),
vecName_(NULL),
timeId_(0),
twoD_(false)
{
 MPI_Comm_rank(comm_,&comm_me_);
 MPI_Comm_size(comm_,&nprocs_);
 
 if(comm_me_==0)
 {
  std::cout << "\n"
  << "***************************************" << "\n"
  << "***              C3PO               ***" << "\n"
  << "*** (c) DCS GmbH, Linz, Austria     ***" << "\n"
  << "*** (c) TU Graz, Graz, Austria      ***" << "\n"
  << "***************************************" << "\n";
 }
  
 MPI_Barrier(comm_);
   
 myC3po_ = new c3po(0,NULL,comm_);
 
 mesh_   = new CSVmesh(myC3po_);
 
 lagrangian_ = new CSVlagrangian(myC3po_,mesh_);
 
 outDirGenerated_ = false;
 
 fO_ = new CSVfieldOperations( mesh_ , fields_ );
 
}

/* ----------------------------------------------------------------------
   c3poCSVInterface Destructors
------------------------------------------------------------------------- */
c3poCSVinterface::~c3poCSVinterface()
{
 delete myC3po_;
 delete mesh_;
}

/* ---------------------------------------------------------------------- */
void c3poCSVinterface::resetAllFields() const
{
 deleteC3POfields();
 delete vecName_;
 
  for(int i=0;i<NofFields_;i++)
     delete csv_[i];
 delete csv_;
 names_.clear();
 Ffieldnames_.clear();
}

/* ---------------------------------------------------------------------- */
void c3poCSVinterface::checkMesh() const
{
 mesh_->checkMesh();
}

/* ---------------------------------------------------------------------- */
void c3poCSVinterface::clearMesh() const
{
 mesh_->clearMesh();
}

/* ---------------------------------------------------------------------- */
void c3poCSVinterface::checkParticles() const
{
 lagrangian_->checkLagrangian();
}

/* ---------------------------------------------------------------------- */
void c3poCSVinterface::clearParticles() const
{
 lagrangian_->deleteParticles();
}

/* ---------------------------------------------------------------------- */
void c3poCSVinterface::createFileList() const
{

 if(comm_me_==0) system("ls data >fileList.c3po");
 
 MPI_Barrier(MPI_COMM_WORLD);
 
 std::ifstream ifnfile ("fileList.c3po",std::ifstream::in);
 
 std::istream& infile_(ifnfile);

 int haveData = 0; 

 while (!infile_.eof())
 { 
  char buf[MAX_CHARS_PER_LINE];
  infile_.getline(buf, MAX_CHARS_PER_LINE);
  
  std::string tmp_("./data/");
  unsigned int sizetmp_=tmp_.size();
  
  tmp_.append(buf);
  
  if(tmp_.size()!=sizetmp_) fileList_.push_back(tmp_);
 
  haveData++;

 }

 if(haveData<2)
 {
     std::cout << "\n proc: " << comm_me_ 
               << " -> Error: Could not read data from directory 'data'. Perhaps it is empty or not in place. \n";

     MPI_Finalize();
     exit(1);
 }

}


/* ---------------------------------------------------------------------- */
void c3poCSVinterface::readInput() const
{ 
 
 char buf[10];
 sprintf(buf,"%i",timeId_);
 
 myC3po_->setTime(buf);
 
 fileName_.assign(fileList_[timeId_]);
 std::ifstream infile_ (fileName_.c_str(),std::ifstream::in);
    
 if (!infile_.is_open())
 {
  std::cout << "\n proc: " << comm_me_ 
            << " -> Error: the file named " 
            << fileName_ << " cannot be opened!\n";
  MPI_Finalize();
  exit(1);
 }
 
 //get field information from CPPPO
 NofVectors_ = myC3po_->getVFnamesNumber();
 
 vecName_= new char*[NofVectors_];
 for (int i=0;i<NofVectors_;i++)
 {
  vecName_[i] = new char[10];
  vecName_[i] = const_cast<char*>(myC3po_->getVFnames(i).c_str());
 }
 
 NofFields_ = NofVectors_*3 + myC3po_->getSFnamesNumber();
 
}


/* ---------------------------------------------------------------------- */
void c3poCSVinterface::parseFile() const
{
  
 std::ifstream infile_ (fileName_.c_str(),std::ifstream::in);
    
 if (!infile_.is_open())
 {
  std::cout << "\n proc: " << comm_me_ << " -> Error: the file named " << fileName_ << " cannot be opened!\n";
  MPI_Finalize();
  exit(1);
 }

 std::istream& infile(infile_);    

 int n=0;
 int s=0;
 std::vector<int> map;


 
   //Parsing the first line, i.e. the names
 
   std::string buf;
   std::vector<std::string> temp_;
   std::getline(infile,buf);
   std::string::iterator it=buf.begin();
   while(it<buf.end())
   {
     std::string str;
     
    while(it!=buf.end())
    {
     if(*it==',') break;
     if(*it==' ')
      {
       it++;
       continue;
      }
     str.push_back(*it);
     it++;    
    }
    temp_.push_back(str); 
    it++;    
   }
   int tsize_=temp_.size();
   bool vectory_ = false;
   bool vectorx_ = false;
   for( int i=0;i<tsize_; i++)
   { 
    for(int j=0;j<NofVectors_;j++)
    {
     if(myC3po_->getVFnames(j).compare(temp_[i])==0)
     {
      
      if(vectorx_)
      {
       std::cout << "\n proc: " << comm_me_ << " -> Error: one dimensional case not supported. Please provide at least two components for every vector field\n";
       MPI_Finalize();
       exit(1);
      }
     
      map.push_back(i);
      names_.push_back(temp_[i]);
      vectorx_=true;
      if(vectory_)
      { 
       twoD_=true;
       std::cout << "\nTwo dimensional case detected! \n";
       vectory_=false;
      } 
         
     }
     else
     {
      unsigned int strlength = temp_[i].size();
     
      if(strlength == myC3po_->getVFnames(j).size())
      {

       for(unsigned int k=0;k<strlength+1;k++)
       {
        if(myC3po_->getVFnames(j)[k] == 'x' &&   temp_[i][k] == 'y'  )
        {
         if(myC3po_->getVFnames(j).compare(temp_[i-1])==0)
         {
          map.push_back(i);
          names_.push_back(temp_[i]);
          vectory_=true;
          vectorx_=false;
         }
         else
         {
          std::cout << "\n proc: " << comm_me_ << " -> Error: Vector field component " << myC3po_->getVFnames(j)<< " and vector field component " << temp_[i] << " are not adjacent in the csv file!\n";
          MPI_Finalize();
          exit(1);
         }
        }
        else if(myC3po_->getVFnames(j)[k] == 'x' &&   temp_[i][k] == 'z'  )
        {

         if(myC3po_->getVFnames(j).compare(temp_[i-2])==0)
         {
          map.push_back(i);
          names_.push_back(temp_[i]);
          vectory_=false;
         }
         else
         {
          std::cout << "\n proc: " << comm_me_ << " -> Error: Vector field component " << myC3po_->getVFnames(j)<< " and vector field component " << temp_[i] << " are not adjacent in the csv file!\n";
          MPI_Finalize();
          exit(1);
         }
        }
        else if(myC3po_->getVFnames(j)[k] != temp_[i][k]   ) break;
       }
      }
     }
    }
    for(int j=0;j<myC3po_->getSFnamesNumber();j++)
    {
     if(myC3po_->getSFnames(j).compare(temp_[i])==0)
     {
      map.push_back(i);
      names_.push_back(temp_[i]);
     
      if(vectorx_)
      {
       std::cout << "\n proc: " << comm_me_ << " -> Error: one dimensional case not supported. Please provide at least two components for every vector field\n";
       MPI_Finalize();
       exit(1);
      }
     
      if(vectory_)
      { 
       twoD_=true;
       std::cout << "\nTwo dimensional case detected! \n";
       vectory_=false;
      } 
          
     }
    }
   }
   
  
   if(vectory_  && !twoD_)
   {
    twoD_=true;
    std::cout << "\nTwo dimensional case detected! \n";
   }

   
   if(vectorx_ )
   {
    std::cout << "\n proc: " << comm_me_ << " -> Error: one dimensional case not supported. Please provide at least two components for every vector field\n";
    MPI_Finalize();
    exit(1);
   } 
    
                 
   int nsize_=names_.size();
   
   if(twoD_)
    if( nsize_!= NofVectors_*2 + myC3po_->getSFnamesNumber() )
    {
     std::cout << "\n proc: " << comm_me_ << " -> Error: Some of the fields you specified in the c3po.json file could not be found.\n";
     MPI_Finalize();
    exit(1);
    }
    
   if(!twoD_)
    if( nsize_!= NofFields_)
    {
     std::cout << "\n proc: " << comm_me_ << " -> Error: Some of the fields you specified in the c3po.json file could not be found.\n";
     MPI_Finalize();
     exit(1);
    }
    
   if(twoD_)
   {
     NofFields_= NofVectors_*2 + myC3po_->getSFnamesNumber();
     csv_ = new double*[NofFields_];
     for(int i=0;i<NofFields_;i++)
      csv_[i] = new double[*(mesh_->NofCells())];
   }
   else
   {
     csv_ = new double*[NofFields_];
     for(int i=0;i<NofFields_;i++)
      csv_[i] = new double[*(mesh_->NofCells())];
   }
        //Skipping those lines that do not belong to this processor 

 while (!infile.eof())
 {

   bool belongs_ = false;
   for(int i=0; i< *(mesh_->NofCells()); i++)
    if(n == *(mesh_->getCellId(i)))
    {
     belongs_=true;
     break;
    }
   if(belongs_) //Storing field values in the csv_ double array
   {
    std::vector<double> temp;
    std::string buf;
    std::getline(infile,buf);
    std::string::iterator it=buf.begin();
    while(it<buf.end())
    {
     std::string str;
     while(it!=buf.end())
     {
      if(*it==',') break;
      str.push_back(*it);
      it++;
     }
     
     temp.push_back(std::atof(str.c_str()));  
     it++;
    }
    
    if(temp.size()==0 && s<*(mesh_->NofCells())-1)
    {
     std::cout << "\n proc: " << comm_me_ << " -> Error: Number of cells in mesh.csv differs from the number of cells in the csv input file!.\n";
     MPI_Finalize();
     exit(1);
    }
    
    if(temp.size()==0 && s==*(mesh_->NofCells())-1)
    {
     break; 
    }
    
    for(int i=0;i<NofFields_;i++)
    {
     int ind_=map[i];
 
     csv_[i][s]=temp[ind_];
    }
    s++;
   }
   else
   {
    std::string buf;
    std::getline(infile,buf);
   } 
  
  n++; 
 } 
 
/* for(int i=0;i<NofFields_;i++)
  std::cout << names_[i] << ",";  
 
 for(int j=0;j<cell_per_proc;j++)
 {
  for(int i=0;i<NofFields_;i++)
   std::cout << csv_[i][j] << ",";
 std::cout << std::endl;
 } */
}


/* ---------------------------------------------------------------------- */
void c3poCSVinterface::registerC3poFields() const
{
 unsigned int numVF = myC3po_->getVFnamesNumber();
 unsigned int numSF = myC3po_->getSFnamesNumber();

 //registering fields in c3po
 for(unsigned int j=0;j<numVF; j++)
  for(int i=0;i<NofFields_;i++)
  {
   if (names_[i].compare(myC3po_->getVFnames(j))==0)
   {
    if(!twoD_)
      myC3po_->GlobalVF(myC3po_->getVFnames(j),csv_[i],csv_[i+1],csv_[i+2]);  
    else
    {
     double * tempz = new double[*(mesh_->NofCells())];
      for(int k=0;k<*(mesh_->NofCells());k++)
       tempz[k]=0.0;
      
     myC3po_->GlobalVF(myC3po_->getVFnames(j),csv_[i],csv_[i+1],tempz);
     dummyZ.push_back(tempz); 
    
    }
   }  
  }
 
 for(unsigned int j=0;j<numSF; j++)
  for(int i=0;i<NofFields_;i++)
  {
   if (names_[i].compare(myC3po_->getSFnames(j))==0)
   {
    myC3po_->GlobalSF(myC3po_->getSFnames(j),csv_[i]);  
   
   }  
  }
}

/* ---------------------------------------------------------------------- */
void c3poCSVinterface::createGradients(int id) const
{
  int OpFilterNum_ = myC3po_->getOpFiltNum();
  std::vector<std::string> GradScalList_ = myC3po_->getGradientScalarList();
 
  int cell_per_proc = *(mesh_->NofCells());  //Get number of cells
  int index = 0;

  int cells [3];
  double cell_size [3];
  double * baseField_ = new double[cell_per_proc];

  // Get number of cells and cell sizes in each coordinate
  for (int c = 0; c<3; c++)
  {
	  cells[c] = (myC3po_->getMaxDomainfromJson(c)-myC3po_->getMinDomainfromJson(c))/myC3po_->cellSizefromJson(c);
	  cell_size[c] = myC3po_->cellSizefromJson(c);
  }
	
  //Create Gradients of Scalar Fields
   for(int n=0;n<OpFilterNum_;n++) //Loop over filtering operations
   {
	for(unsigned int sf=0; sf<(myC3po_->scalarFTF(n)).size();sf++)  //Loop over scalar fields for filtering in each operation
	{
		for(unsigned int scal = 0; scal < GradScalList_.size(); scal++)  // Loop over the scalar fields specified for gradient calculations
		{
		std::string name("grad");

	    //look for the field name
		if(GradScalList_[scal].compare((myC3po_->scalarFTF(n))[sf])==0)
		{	
			int check = check_repetition(n, sf);  //Check if the gradients has already been calculated for the current field

			if (check == 0)
			{
				name.append(GradScalList_[scal].c_str());

				double * gradx = new double[cell_per_proc];
				double * grady = new double[cell_per_proc];
				double * gradz = new double[cell_per_proc];
				
				std::string current(GradScalList_[scal].c_str());
				current.append("_");
				current.append(myC3po_->getOpFilterName(n)); 
				current.append("_");
				current.append((myC3po_->getFilterName(id)).c_str()); 
				for (unsigned int j = 0; j < Ffieldnames_.size(); j++ )  // Get the correct ID for the field for gradient calculations
				{					
					if (Ffieldnames_[j].compare(current)==0)
					{
						index = j - id*NofFiltFields_;
						break;
					}
				}

				for (int idC=0; idC<*(mesh_->NofCells()); idC++)  // Assign values for the field of interest
				{
						* (baseField_+idC) = fields_[index][idC];
				}

				double delx[cell_per_proc], dely[cell_per_proc], delz[cell_per_proc];
		  
				fO_->evaluateGradient(sf,gradx,grady,gradz, cells, cell_size, baseField_);  // Calculate gradients

				for (int idC=0; idC<*(mesh_->NofCells()); idC++)  //Save gradient values
				{
						delx[idC] = * (gradx+idC); 
						dely[idC] = * (grady+idC); 
						delz[idC] = * (gradz+idC); 
				}

				myC3po_->registerVF(name.c_str(),gradx,grady,gradz);

				for (int idC=0; idC<*(mesh_->NofCells()); idC++)  //Reassign gradient values
				{
						* (gradx+idC) = delx[idC];
						* (grady+idC) = dely[idC];
						* (gradz+idC) = delz[idC];
				}
			 
				gradientFields_.push_back(gradx);
				gradientFields_.push_back(grady);
				gradientFields_.push_back(gradz);
			}			
		}		
	}
   } 
  }
}

/* ---------------------------------------------------------------------- */
void c3poCSVinterface::printCSV(int id) const
{
 
 if(!myC3po_->writeFields()) return;
 
 std::string filterName_(myC3po_->getFilterName(id));
 std::string fileName_("./c3po_CSV_fields/");
 
 //create directory
 if( (!outDirGenerated_) && (comm_me_==0) )
 {
        if(stat(fileName_.c_str(), &dirLogicCSV) == 0 && S_ISDIR(dirLogicCSV.st_mode))
            std::cout << ""; 
        else
        {
            std::string command("mkdir "+fileName_);
            system(command.c_str());
        }
        outDirGenerated_ = true;
 }
 MPI_Barrier(MPI_COMM_WORLD); //must sync here!
 
 
 fileName_.append(filterName_);
 char buf[80];
 sprintf(buf,"_time%i.csv",timeId_);
 fileName_.append(buf);
 
 if(comm_me_==0)
 {
   std::ofstream outfile_(fileName_.c_str());
  
   //writing the header
   for(int i=0;i<NofFiltFields_+NofFiltVarianceFields_;i++)
   {
    if(i==NofFiltFields_+NofFiltVarianceFields_-1)
     outfile_ << Ffieldnames_[i+id*(NofFiltFields_+NofFiltVarianceFields_)] 
              << std::endl;
    else
     outfile_ << Ffieldnames_[i+id*(NofFiltFields_+NofFiltVarianceFields_)]  
              << ",";
   }
 }
 
 MPI_Barrier(MPI_COMM_WORLD);
 
 int n=0;
 while(n<*(mesh_->NofCellsGlobal()))
 {
  bool belongs_ = false;
  int id_;
  for(int i=0; i< *(mesh_->NofCells()); i++)
   if(n == *(mesh_->getCellId(i)))
   {
    belongs_=true;
    id_=i;
    break;
   }
  
  if(belongs_)
  {
   std::ofstream outfile_(fileName_.c_str(), std::ios::out | std::ios::app);
  
   //writing data
   for(int i=0;i<NofFiltFields_+NofFiltVarianceFields_;i++)
   {
    if(i==NofFiltFields_+NofFiltVarianceFields_-1)
     outfile_ << fields_[i][id_] << std::endl;
    else
     outfile_ << fields_[i][id_] << ",";
   }
   outfile_.close();
  }

  //other processors must wait, the file is written sequentially
  MPI_Barrier(MPI_COMM_WORLD);
  n++;
 }

}


/* ---------------------------------------------------------------------- */
void c3poCSVinterface::deleteFields() const
{
 myC3po_->resetFields();
 
 for (int i=0;i<(NofFiltFields_+NofFiltVarianceFields_);i++)
  delete fields_[i];
 delete fields_;
 
 for(unsigned int i=0; i<gradientFields_.size(); i++)
  delete gradientFields_[i];
  
 gradientFields_.clear();
}


/* ---------------------------------------------------------------------- */
void c3poCSVinterface::createFilterFields(int id) const
{
 std::string filterName_(myC3po_->getFilterName(id));
 int OpFilterNum_ = myC3po_->getOpFiltNum();
 NofFiltFields_=0;
 NofFiltVarianceFields_=0;
 int cell_per_proc = *(mesh_->NofCells());
 
 for(int i=0;i<OpFilterNum_;i++)
 {
  NofFiltFields_ += 3*myC3po_->vectorFTF(i).size() 
                  +   myC3po_->scalarFTF(i).size();
  NofFiltVarianceFields_ += 3*myC3po_->vectorFTFVariance(i).size() 
                         +    myC3po_->scalarFTFVariance(i).size();
 }
 
 //creating arrays for filtering operations 
 fields_ = new double*[NofFiltFields_+NofFiltVarianceFields_];
 for (int i=0;i<(NofFiltFields_+NofFiltVarianceFields_);i++)
  fields_[i] = new double[cell_per_proc];
  
 int s=0;
 for(int n=0;n<OpFilterNum_;n++)
 {

  std::string OpName_(myC3po_->getOpFilterName(n));

  //Filtered fields
  for( int j=0;j<NofFields_; j++)
   for(unsigned int i=0;i<(myC3po_->vectorFTF(n)).size();i++)
    if (names_[j].compare((myC3po_->vectorFTF(n))[i])==0)
    {
      std::string fieldName_(names_[j]); 
      fieldName_.append("_");
      fieldName_.append(OpName_.c_str()); 
      fieldName_.append("_");
      fieldName_.append(filterName_.c_str()); 
      Ffieldnames_.push_back(fieldName_+"_x");
      Ffieldnames_.push_back(fieldName_+"_y");
      Ffieldnames_.push_back(fieldName_+"_z"); 
      myC3po_->registerVF(names_[j],fields_[s],fields_[s+1],fields_[s+2]);
      s+=3;
    }
  
 for( int j=0;j<NofFields_; j++)
   for(unsigned int i=0;i<(myC3po_->scalarFTF(n)).size();i++)
    if (names_[j].compare((myC3po_->scalarFTF(n))[i])==0)
    {
      std::string fieldName_(names_[j]); 
      fieldName_.append("_");
      fieldName_.append(OpName_.c_str()); 
      fieldName_.append("_");
      fieldName_.append(filterName_.c_str()); 
      Ffieldnames_.push_back(fieldName_);      
      myC3po_->registerSF(names_[j],fields_[s]);
      s++;
    }
    
  //Variance fields
  for( int j=0;j<NofFields_; j++)
   for(unsigned int i=0;i<(myC3po_->vectorFTFVariance(n)).size();i++)
    if (names_[j].compare((myC3po_->vectorFTFVariance(n))[i])==0)
    {
      std::string fieldName_(names_[j]); 
      fieldName_.append("_");
      fieldName_.append(OpName_.c_str()); 
      fieldName_.append("_variance");
      std::ostringstream Str; Str << i;  fieldName_.append(Str.str());
      fieldName_.append("_");
      fieldName_.append(filterName_.c_str()); 
      Ffieldnames_.push_back(fieldName_+"_x");
      Ffieldnames_.push_back(fieldName_+"_y");
      Ffieldnames_.push_back(fieldName_+"_z"); 
      myC3po_->registerVF(names_[j],fields_[s],fields_[s+1],fields_[s+2]);
      s+=3;
    }
    
  for( int j=0;j<NofFields_; j++)
   for(unsigned int i=0;i<(myC3po_->scalarFTFVariance(n)).size();i++)
    if (names_[j].compare((myC3po_->scalarFTFVariance(n))[i])==0)
    {
      std::string fieldName_(names_[j]); 
      fieldName_.append("_");
      fieldName_.append(OpName_.c_str()); 
      fieldName_.append("_variance");
      std::ostringstream Str; Str << i;  fieldName_.append(Str.str());
      fieldName_.append("_");
      fieldName_.append(filterName_.c_str()); 
      Ffieldnames_.push_back(fieldName_);      
      myC3po_->registerSF(names_[j],fields_[s]);
      s++;
    }

 
 }
}


/* ---------------------------------------------------------------------- */
void c3poCSVinterface::runFilter(int id) const
{
 
  createFilterFields(id);    
  myC3po_->runFilters(id);
  createGradients(id);
     
}


/* ---------------------------------------------------------------------- */
void c3poCSVinterface::runSampling() const
{
  myC3po_->runSampling();
}


/* ---------------------------------------------------------------------- */
void c3poCSVinterface::runBinning() const
{
  myC3po_->runBinning();
}


/* ---------------------------------------------------------------------- */
void c3poCSVinterface::runC3po() const
{ 
 myC3po_->preRunOperations();
 
 int NofFilters_=myC3po_->numberOfFilters();
 for(int id=0;id<NofFilters_;id++)
  {
   runFilter(id);
   printCSV(id);
   runSampling();
   runBinning();
   deleteFields(); 
  }
 myC3po_->postRunOperations();
}


/* ---------------------------------------------------------------------- */
void c3poCSVinterface::deleteC3POfields() const
{
 myC3po_->resetGlobalFields();
 for(unsigned int i=0;i<dummyZ.size();i++)
  delete dummyZ[i];
  
 dummyZ.clear();
} 


/* ---------------------------------------------------------------------- */
void c3poCSVinterface::FLUENTrunC3PO() const
{

 createFileList();
 for(unsigned int time=0;time<fileList_.size();time++) 
 {
 
  timeId_=time;
  
  if(comm_me_==0)
   printf("\n===> CPPPO is processing time %i/%lu from file %s  \n",
           timeId_,fileList_.size() -1,fileList_[time].c_str());
  
  readInput();
  parseFile();
  registerC3poFields();
  runC3po();
  resetAllFields();

 }
}

 /* ---------------------------------------------------------------------- */
 int c3poCSVinterface::check_repetition (int n, int sf) const
 {
	 		int check = 0;

			if (n == 0)
			{
				for (int sf_test=0; sf_test<sf; sf_test++)
				{
					if ((myC3po_->scalarFTF(0))[sf_test].compare((myC3po_->scalarFTF(n))[sf])==0)
					{
						check = 1;	
					}				
				}
			}
			int n_test = 0 ;
			for (n_test = 0; n_test < n; n_test++)
			{
				for (int sf_test=0; sf_test<(myC3po_->scalarFTF(n_test)).size(); sf_test++)
				{
					if ((myC3po_->scalarFTF(n_test))[sf_test].compare((myC3po_->scalarFTF(n))[sf])==0)
					{
						check = 1;
					}
				}
			}

			if (n_test == n && n_test != 0)
			{
				for (int sf_test=0; sf_test<sf; sf_test++)
				{
					if ((myC3po_->scalarFTF(n_test))[sf_test].compare((myC3po_->scalarFTF(n))[sf])==0)
					{
						check = 1;	
					}				
				}
			}

			return check; //check = 1 if a repeated filtered field is found
 }


