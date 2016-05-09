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
#include "region.h"
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

region::region(std::string regionName, std::vector<std::string> filterNames, bool readCellsFromJson, c3po *ptr)
:
c3poBaseAccessible(ptr),
regionName_(regionName),
regionId_(-1),
filterNames_(filterNames),
readCellsFromJson_(readCellsFromJson),
volume_(0),
surfaceArea_(0)
{


}

//--------------------------------------------------------------------------------------
region::~region()
{
 delete nParticlesProc_;
}

//--------------------------------------------------------------------------------------

/* ----------------------------------------------------------------------
   Add  cells to memory
------------------------------------------------------------------------- */
void region::addCell(int & IDInternal, bool isOnSurface)
{
 
 if(isOnSurface)
  cellIDsSurface_.push_back(IDInternal);
 else
  cellIDsInternal_.push_back(IDInternal);
 
}
/* ----------------------------------------------------------------------
   Clear region content
------------------------------------------------------------------------- */
void region::clearCells()
{
  cellIDsSurface_.clear();
  cellIDsInternal_.clear();
}
/* ----------------------------------------------------------------------
   Calculate total volume of region
------------------------------------------------------------------------- */
void region::evaluateVolume()
{
 volume_=0;

 //evaluate internal volume
 for(unsigned int cell=0;cell<cellIDsInternal_.size();cell++)
  volume_ += *mesh().CellVol( cellIDsInternal_[cell] ); 

}
/* ----------------------------------------------------------------------
   Calculate total surface area of region (this is not a real area!)
------------------------------------------------------------------------- */
void region::evaluateSurfaceArea()
{
 
 surfaceArea_=0;
 
 //evaluate surface volume
 for(unsigned int cell=0;cell<cellIDsSurface_.size();cell++)
  surfaceArea_ += *mesh().CellVol( cellIDsSurface_[cell] ); 

 
}
/* ----------------------------------------------------------------------
   Calculate center of mass (constant density!)
------------------------------------------------------------------------- */
void region::evaluateCenterOfMass()
{
 
 evaluateVolume();
 
 centerOfMass_[0]=0;
 centerOfMass_[1]=0;
 centerOfMass_[2]=0;
  
 for(unsigned int cell=0;cell<cellIDsInternal_.size();cell++)
 {
  centerOfMass_[0] += *mesh().CellVol( cellIDsInternal_[cell] ) * ( mesh().CellCentre( cellIDsInternal_[cell] )[0]) ; 
  centerOfMass_[1] += *mesh().CellVol( cellIDsInternal_[cell] ) * ( mesh().CellCentre( cellIDsInternal_[cell] )[1]) ; 
  centerOfMass_[2] += *mesh().CellVol( cellIDsInternal_[cell] ) * ( mesh().CellCentre( cellIDsInternal_[cell] )[2]) ; 
 }

 centerOfMass_[0] /= volume_;
 centerOfMass_[1] /= volume_;
 centerOfMass_[2] /= volume_;
 
}
/* ----------------------------------------------------------------------
   Write region to disk
------------------------------------------------------------------------- */
void region::write()
{
 std::string filename_("./c3po_dataStorage_regions/");
 filename_.append(regionName_);
 
 evaluateVolume();
 evaluateSurfaceArea();
 evaluateCenterOfMass();
 
 #ifdef H5_LIB
 
 if(input().dumpFormat().compare("hdf5")==0)
 {
  std::string file(filename_);
  
  file.append("_cellList");
  file.append(".h5");
  
  OneArrayToH5(file, "internalCells"       , &cellIDsInternal_[0], cellIDsInternal_.size());
  OneArrayToH5(file, "surfaceCells"        , &cellIDsSurface_[0] , cellIDsSurface_.size() );
  
  file.assign(filename_);
  file.append(".h5");
  OneArrayToH5(file, "regionVolume"        , &volume_            , 1                      );
  OneArrayToH5(file, "surfaceArea"         , &surfaceArea_       , 1                      );
  OneArrayToH5(file, "centerOfMass (x,y,z)", &centerOfMass_[0]    , 3                      );
  
 
 }
 
 #endif
 
 if(input().dumpFormat().compare("json")==0)
 {
  std::string file(filename_);
  
  file.append("_cellList");
  file.append(".json");
  
  std::vector< int* >           dataVec_;
  std::vector<std::string >    dataName_;
  std::vector<int>              dataNum_;
  
  dataVec_.push_back( &cellIDsInternal_[0] );
  dataName_.push_back( "internalCells" );
  dataNum_.push_back(cellIDsInternal_.size());
 
  dataVec_.push_back( &cellIDsSurface_[0] );
  dataName_.push_back( "surfaceCells" );
  dataNum_.push_back(cellIDsSurface_.size());
  
  output().createQJsonArrays(file,"cellList",dataName_,dataVec_, -1, true, &dataNum_); 
 
  std::vector< double* >       dataVecD_;
  dataName_.clear();
  dataNum_.clear();
  
  file.assign(filename_);
  file.append(".json");
   
  dataVecD_.push_back( &volume_ );
  dataName_.push_back( "regionVolume" );
  dataNum_.push_back(1);
 
  dataVecD_.push_back( &surfaceArea_ );
  dataName_.push_back( "surfaceArea" );
  dataNum_.push_back(1);
  
  dataVecD_.push_back( &centerOfMass_[0] );
  dataName_.push_back( "centerOfMass" );
  dataNum_.push_back(3);
  
  output().createQJsonArrays(file,"regionData",dataName_,dataVecD_, -1, true, &dataNum_); 
 }
 
}
///* ---------------------------------------------------------------------------*/
//void region::readParticles() const
//{
// 
// timer().stamp();

// if(!readParticlesFromJson_) return;
// 
// readParticlesFromJson_=false;
// 
// input().readParticles(&parPos_, groupName_);
// 
// for(unsigned int par=0;par<parPos_.size()/3;par++)
// {
// 
//  //Check if within processor boundaries
//  if( ! (
//         ( parPos_[3*par] < mesh().dLocalMax()[0]  )  
//          &&
//          
//         ( parPos_[3*par] > mesh().dLocalMin()[0]  ) 
//        )
//    ) continue;
// 
//  if( ! (
//         ( parPos_[3*par+1] < mesh().dLocalMax()[1]  )  
//          &&
//          
//         ( parPos_[3*par+1] > mesh().dLocalMin()[1]  ) 
//        )
//    ) continue;

//  if( ! (
//         ( parPos_[3*par+2] < mesh().dLocalMax()[2]  )  
//          &&
//          
//         ( parPos_[3*par+2] > mesh().dLocalMin()[2]  ) 
//        )
//    ) continue;

// 
//  Particle* par_ = new Particle();
//  particles_.push_back(par_);
//  
//  par_->setpos(&parPos_[3*par]);
//  par_->setId(par);
// 
// }
// 
//// removeGhostParticles(); This function will badly affect parallel scalability!!
// timer().stamp(TIME_INPUT);
//}

