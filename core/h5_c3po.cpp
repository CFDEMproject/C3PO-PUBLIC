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
#ifdef H5_LIB

#include "psctype.h"
#include "stdlib.h"
#include "hdf5.h"
#include "h5_c3po.h"
#include "H5Cpp.h"
#include <string>

#define RANK 2

using namespace H5;
using namespace C3PO_NS;

using namespace H5_C3PO_NS;

// * * * * * * * * * * * * * * * * * * * * * * * * * *
void H5_C3PO_NS::createH5file(std::string filename_) 
{
    H5File* file_ = new H5File( filename_.c_str(), H5F_ACC_TRUNC );
    file_->close();
    delete file_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * *
void H5_C3PO_NS::createExtendibleDataset(std::string FILE_NAME,const char* datasetName_)
{

   hsize_t dims[2] = { 0, 1}; // dataset dimensions at creation
   hsize_t maxdims[2] = {H5S_UNLIMITED, H5S_UNLIMITED};
   DataSpace mspace1( RANK, dims, maxdims);

  H5File* file=new H5File( FILE_NAME.c_str(),H5F_ACC_RDWR );
  IntType datatype( PredType::NATIVE_DOUBLE );         //Define datatype for the data
  datatype.setOrder( H5T_ORDER_LE );
  
  DSetCreatPropList cparms;
  hsize_t chunk_dims[2] ={6, 1};
  cparms.setChunk( RANK, chunk_dims );
  
  //Set fill value for the dataset
  
  int fill_val = 1.0;
  cparms.setFillValue( PredType::NATIVE_DOUBLE, &fill_val);
  
  DataSet dataset = file->createDataSet( datasetName_, PredType::NATIVE_DOUBLE, mspace1, cparms);
  
  file->close();
  
  delete file;

}

// * * * * * * * * * * * * * * * * * * * * * * * * * *
void H5_C3PO_NS::addOneArrayToH5(std::string FILE_NAME, const char* datasetName_, double* data, int NX)
{
  

  hid_t file, dset,dspace;
  file = H5Fopen (FILE_NAME.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  dset = H5Dopen (file, datasetName_, H5P_DEFAULT);
  dspace = H5Dget_space(dset);
  hsize_t dims[RANK];
  H5Sget_simple_extent_dims(dspace, dims, NULL);
  
  H5Dclose(dset);
  H5Fclose(file);
  
  H5File* file_=new H5File( FILE_NAME.c_str(),H5F_ACC_RDWR );
  IntType datatype( PredType::NATIVE_DOUBLE );         //Define datatype for the data
  datatype.setOrder( H5T_ORDER_LE );
  
  DataSet* dataset = new DataSet(file_->openDataSet( datasetName_));
  
  hsize_t ext[2]={hsize_t(NX),hsize_t(0)};
  hsize_t newdims[2]={dims[0]+ext[0],1};
  dataset->extend( newdims );
  
  DataSpace fspace3 = dataset->getSpace();
  fspace3.selectHyperslab( H5S_SELECT_SET, ext, dims );

   /*
  * Define memory space.
  */
  DataSpace mspace3( RANK, ext );
  /*
  * Write the data to the hyperslab.
  */
  dataset->write( data, PredType::NATIVE_DOUBLE, mspace3, fspace3 );
  
  delete dataset;   
  file_->close();
  delete file_; 

}

// * * * * * * * * * * * * * * * * * * * * * * * * * *
void H5_C3PO_NS::addOneArrayToH5(std::string FILE_NAME, const char* datasetName_, int* data, int NX)
{

  hid_t file, dset,dspace;
  file = H5Fopen (FILE_NAME.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  dset = H5Dopen (file, datasetName_, H5P_DEFAULT);
  dspace = H5Dget_space(dset);
  hsize_t dims[RANK];
  H5Sget_simple_extent_dims(dspace, dims, NULL);
  
  H5Dclose(dset);
  H5Fclose(file);
  
  H5File* file_=new H5File( FILE_NAME.c_str(),H5F_ACC_RDWR );
  IntType datatype( PredType::NATIVE_INT );         //Define datatype for the data
  datatype.setOrder( H5T_ORDER_LE );
  
  DataSet* dataset = new DataSet(file_->openDataSet( datasetName_));
  
  hsize_t ext[2]={hsize_t(NX),hsize_t(0)};
  hsize_t newdims[2]={dims[0]+ext[0],1};
  dataset->extend( newdims );
  
  DataSpace fspace3 = dataset->getSpace();
  fspace3.selectHyperslab( H5S_SELECT_SET, ext, dims );

   /*
  * Define memory space.
  */
  DataSpace mspace3( RANK, ext );
  /*
  * Write the data to the hyperslab.
  */
  dataset->write( data, PredType::NATIVE_INT, mspace3, fspace3 );
  
  delete dataset;   
  file_->close();
  delete file_; 

}

// * * * * * * * * * * * * * * * * * * * * * * * * * *
void H5_C3PO_NS::OneArrayToH5(std::string FILE_NAME, const char* datasetName_, double* data, int NX) 
{
  

  H5File* file=new H5File( FILE_NAME.c_str(),H5F_ACC_RDWR );
  IntType datatype( PredType::NATIVE_DOUBLE );         //Define datatype for the data
  datatype.setOrder( H5T_ORDER_LE );
  
  const int NY=1;
  hsize_t     dimsf[2]={hsize_t(NX),hsize_t(NY)};              /* dataset dimensions */      
  DataSpace dataspace( 2, dimsf);
  
  DataSet* dataset = new DataSet(file->createDataSet(datasetName_, datatype, dataspace )); //Create a new dataset within the file using defined dataspace
 
  dataset->write( data, PredType::NATIVE_DOUBLE ); //Write the data to the dataset using default memory space, file

  delete dataset;   
  file->close();
  delete file; 


}

// * * * * * * * * * * * * * * * * * * * * * * * * * *
void H5_C3PO_NS::OneArrayToH5(std::string FILE_NAME, const char* datasetName_, int* data, int NX)
{

  H5File* file=new H5File( FILE_NAME.c_str(),H5F_ACC_RDWR );
  IntType datatype( PredType::NATIVE_INT );         //Define datatype for the data
  datatype.setOrder( H5T_ORDER_LE );
  
  const int NY=1;
  hsize_t     dimsf[2]={hsize_t(NX),hsize_t(NY)};              /* dataset dimensions */      
  DataSpace dataspace( 2, dimsf);
  
  DataSet* dataset = new DataSet(file->createDataSet(datasetName_, datatype, dataspace )); //Create a new dataset within the file using defined dataspace
 
  dataset->write( data, PredType::NATIVE_INT ); //Write the data to the dataset using default memory space, file

  delete dataset;    
  file->close(); 
   delete file; 
}

// * * * * * * * * * * * * * * * * * * * * * * * * * *
void H5_C3PO_NS::TwoArrayToH5(std::string FILE_NAME, const char* datasetName_, double** dataXY, int NX) 
{

  H5File* file=new H5File( FILE_NAME.c_str(),H5F_ACC_RDWR );
  IntType datatype( PredType::NATIVE_DOUBLE );         //Define datatype for the data
  datatype.setOrder( H5T_ORDER_LE );
  
  const int NY=2;
  hsize_t     dimsf[2]={hsize_t(NX),hsize_t(NY)};              /* dataset dimensions */      
  DataSpace dataspace( 2, dimsf);
  
  DataSet* dataset = new DataSet(file->createDataSet(datasetName_, datatype, dataspace )); //Create a new dataset within the file using defined dataspace
 
  dataset->write( dataXY, PredType::NATIVE_DOUBLE ); //Write the data to the dataset using default memory space, file

  delete dataset;    
  file->close(); 
   delete file; 


}

// * * * * * * * * * * * * * * * * * * * * * * * * * *
void H5_C3PO_NS::TwoArrayToH5(std::string FILE_NAME, const char* datasetName_, double dataXY[][2], int NX) 
{

  H5File* file=new H5File( FILE_NAME.c_str(),H5F_ACC_RDWR );
  IntType datatype( PredType::NATIVE_DOUBLE );         //Define datatype for the data
  datatype.setOrder( H5T_ORDER_LE );
  
  const int NY=2;
  hsize_t     dimsf[2]={hsize_t(NX),hsize_t(NY)};              /* dataset dimensions */      
  DataSpace dataspace( 2, dimsf);
  
  DataSet* dataset = new DataSet(file->createDataSet(datasetName_, datatype, dataspace )); //Create a new dataset within the file using defined dataspace
 
  dataset->write( dataXY, PredType::NATIVE_DOUBLE ); //Write the data to the dataset using default memory space, file

  delete dataset;    
  file->close(); 
   delete file; 


}

// * * * * * * * * * * * * * * * * * * * * * * * * * *
void H5_C3PO_NS::ThreeArrayToH5(std::string FILE_NAME, const char* datasetName_, double dataXY[][3], int NX) 
{

  H5File* file=new H5File( FILE_NAME.c_str(),H5F_ACC_RDWR );
  IntType datatype( PredType::NATIVE_DOUBLE );         //Define datatype for the data
  datatype.setOrder( H5T_ORDER_LE );
  
  const int NY=3;
  hsize_t     dimsf[2]={hsize_t(NX),hsize_t(NY)};              /* dataset dimensions */      
  DataSpace dataspace( 2, dimsf);
  
  DataSet* dataset = new DataSet(file->createDataSet(datasetName_, datatype, dataspace )); //Create a new dataset within the file using defined dataspace
 
  dataset->write( dataXY, PredType::NATIVE_DOUBLE ); //Write the data to the dataset using default memory space, file

  delete dataset;    
  file->close(); 
   delete file; 

}

// * * * * * * * * * * * * * * * * * * * * * * * * * *

#endif

