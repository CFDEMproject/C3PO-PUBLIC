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

#include "mesh_check.H"
#include "fvMesh.H"
#include "polyMesh.H"
#include "c3po.h"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include <cmath>
#include <algorithm> 

using namespace C3PO_NS;
using namespace std;

meshCheck::meshCheck(const Foam::fvMesh& mesh,c3po* ptr)
:
mesh_(mesh),
C3po_(ptr)
{}

meshCheck::~meshCheck() {}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void meshCheck::checkMyMesh() const
{
    
    C3po_->clearMesh();
    
    double cell_num;
    bool struct_;
   
    cell_num = returnReduce(mesh_.nCells(), sumOp<label>());                 //returns the number of cells in mesh_
    
    C3po_->setNofCells(mesh_.nCells());
    
     //For structured mesh only!!
    determine_cellsize();
    
    compute_globals(); 
    registerC3POcells();
	
   
	
    double tolerance = C3po_->meshCheckTolerance();
    if( C3po_->meshVerbose() )
    {
       // MPI_Barrier(communicator_);
        Pout << "selected mesh_ in proc "<< " has " << cell_num << " cells!" << endl;
        Info << "Now checking if cell volume is within " << tolerance << " of tolerance" << endl;
    }

    double vol;
    double vol_ges=0;
    double h=0;
    double vol_avg_time_max, vol_avg_time_min;
    
    forAll(mesh_.cells(), celli )
    {
        h++;                                                              // counter for evaluated cells
        vol=mesh_.V()[celli];                                               // returns actual volume for each celli
        vol_ges+=vol;                                                       // summation of each volume
        vol_avg_time_max=(vol_ges/h)+((vol_ges/h)*tolerance);               // setting maximum and minimum reference cell volume
        vol_avg_time_min=(vol_ges/h)-((vol_ges/h)*tolerance);

        if (C3po_->meshVerbose())
        {
            Info << "current cells volume = " << vol << endl;
            Info<< "volume average maximum = " << vol_avg_time_max << endl;
            Info<< "volume average minumum = " << vol_avg_time_min << endl;
            Info << "current volume = " << vol<<  ", overall volume = "<< vol_ges<< endl;
            Info << "current number of cells looped over = " << h << endl;
        }

        if (vol >= vol_avg_time_min && vol <= vol_avg_time_max)             // true if cell volume is equal
        {
            h=h;                                                            // loop gets the right counter for cells
        }
        else
        {
            h=0;      														// h can never be equal to numbmer of cells in this case
        }
    }

    if (h==(number_of_cells_[0]*number_of_cells_[1]*number_of_cells_[2]))
    {
        Pout << "meshCheck::checkMyMesh: mesh sucsessfully tested - all volumes are equal. You can use IJK selectors for this mesh." << endl;
        struct_=true;
    }
    else
    {
        Info << "meshCheck::checkMyMesh: mesh sucsessfully tested - WARNING: volumes are NOT equal. Do not use IJK selectors for this mesh."<<   endl;
        struct_=false;
	
    } 
    
    C3po_->checkIJK(struct_);
}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
void meshCheck::registerC3POcells() const
{
 
// forAll(mesh_.cells(), celli )
// {
//  C3po_->registerCell(mesh_.V()[celli],&(mesh_.C()[celli].component(0)),&(mesh_.C()[celli].component(1)),&(mesh_.C()[celli].component(2)));
 
// }

   C3po_->registerCells( const_cast<double*>(&(mesh_.V()[0])),const_cast<double*>( &(mesh_.C()[0].component(0)) ),mesh_.nCells(),3) ;
}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void meshCheck::determine_cellsize() const
{
  if( C3po_->meshVerbose() )
	Info << "meshCheck::determine_cellsize: determining cell and mesh size now" << endl;

	//holds the mesh length in each direction [m], only centroids are considered!
	mesh_length_[0] = max(mesh_.cellCentres().component(0)) - min(mesh_.cellCentres().component(0));	//x
	mesh_length_[1] = max(mesh_.cellCentres().component(1)) - min(mesh_.cellCentres().component(1));	//y
	mesh_length_[2] = max(mesh_.cellCentres().component(2)) - min(mesh_.cellCentres().component(2));	//z
	
	//holds the middle point of the first cell - after running checkMyMesh this is be the half cell width in each direction
	//Be careful when using cyclic BCs 
	// cell_middle_cell_0_[0]=abs( (mesh_.C()[0].component(0))- (min(mesh_.points().component(0))));	//x		
        // cell_middle_cell_0_[1]=abs( (mesh_.C()[0].component(1))- (min(mesh_.points().component(1))));	//y
        // cell_middle_cell_0_[2]=abs( (mesh_.C()[0].component(2))- (min(mesh_.points().component(2))));	//z
      
        //holds the length of the first cell in each direction (determined by the first cell)
	// cell_size_[0]= 2.0 * cell_middle_cell_0_[0];	//x		
        // cell_size_[1]= 2.0 * cell_middle_cell_0_[1];	//y
        // cell_size_[2]= 2.0 * cell_middle_cell_0_[2];	//z
        
        cell_size_[0] = 0.;
        cell_size_[1] = 0.;
        cell_size_[2] = 0.;
        
        const faceList & ff = mesh_.faces();
        const pointField & pp = mesh_.points();
        
        const cell & cc = mesh_.cells()[0];
        labelList pLabels(cc.labels(ff));
        
        //Use Cell 0 to get the cell size
        forAll (pLabels, pointi)
        {
         
         point pLocal = pp[pLabels[pointi]];
         cell_size_[0] += abs( pLocal.component(0) - mesh_.C()[0].component(0) )/4;
         cell_size_[1] += abs( pLocal.component(1) - mesh_.C()[0].component(1) )/4;
         cell_size_[2] += abs( pLocal.component(2) - mesh_.C()[0].component(2) )/4;
        
        }
        
	
        mesh_length_[0] += cell_size_[0];
        mesh_length_[1] += cell_size_[1];
        mesh_length_[2] += cell_size_[2];
        
        
	//holds the number of cells in each direction by dividing the length of the mesh in each direction by the cell width of the first cell
	number_of_cells_[0]=lround(mesh_length_[0]/cell_size_[0]);	//x
	number_of_cells_[1]=lround(mesh_length_[1]/cell_size_[1]);	//y  
	number_of_cells_[2]=lround(mesh_length_[2]/cell_size_[2]);	//z
    
  
    
    MPI_Barrier(MPI_COMM_WORLD);
    C3po_->setCells(&number_of_cells_[0],&number_of_cells_global_[0], &cell_size_[0]);
 	
} 

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void meshCheck::compute_globals() const
{ 
    maxDomain[0]=max(mesh_.points().component(0)) ;  
    maxDomain[1]=max(mesh_.points().component(1)) ;
    maxDomain[2]=max(mesh_.points().component(2)) ;
    
    minDomain[0]=min(mesh_.points().component(0)) ;  
    minDomain[1]=min(mesh_.points().component(1)) ;
    minDomain[2]=min(mesh_.points().component(2)) ;
    
    MPI_Allreduce(maxDomain,maxDomainGlobal,3,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(minDomain,minDomainGlobal,3,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    
    C3po_->registerDomainInfo(maxDomain,minDomain,maxDomainGlobal,minDomainGlobal);
}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

