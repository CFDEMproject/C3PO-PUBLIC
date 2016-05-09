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

#include "CSVfieldOperations.h"
#include "mesh.h"

using namespace C3PO_NS;

/* ----------------------------------------------------------------------
   fieldOperations Constructors
------------------------------------------------------------------------- */
CSVfieldOperations::CSVfieldOperations( CSVmesh  *        meshPtr, 
                                        double  **&     fieldsPtr
                                      )
:
mesh_(meshPtr),
initialized_(true)
{
}

/* ----------------------------------------------------------------------
   fieldOperations Destructors
------------------------------------------------------------------------- */
CSVfieldOperations::~CSVfieldOperations()
{
}

/* ----------------------------------------------------------------------
   locating functions (IJK only)
------------------------------------------------------------------------- */
void CSVfieldOperations::cellIDtoIJK(int cellID_, int * IJK, int cells [3] )
{
 
 int x = cells[0];    
 int y = cells[1]/**(cell_global_y())*/;
 IJK[2] = int(cellID_ / (x*y)); //estimation of z-layer
 IJK[1] = int((cellID_-(x*y*IJK[2]))/x); //estimation of y-layer
 IJK[0] = int((cellID_-(x*y*IJK[2]))-((IJK[1]*x))); //estimation of x-layer
}

/*------------------------------------------------------------------------- */
int CSVfieldOperations::IJKtoCellID(int * IJK, int cells [3])
{
 
  int CellID_OF;          
   
    CellID_OF =   (IJK[2])*( cells[0] * cells[1] )
                + (IJK[1])* cells[0]
                +  IJK[0];
    
    return CellID_OF;

}

/*------------------------------------------------------------------------- */
void CSVfieldOperations::evaluateGradient( int      fieldID_ ,
                       double *     gradx,
                       double *     grady, 
                       double *     gradz,
					   int cells [3],
					   double cell_size [3],
					   double * baseField_
                     )
{
 //INPUT: 
 // baseField_    ... id of the field to be sampled (field MUST be scalar)
 //OUTPUT: 
 // gradx       ... field containing the x-components of the gradient (length = length of field)

 int    IJKCurrCell[3];
 int    IJKMinusNeighbor[3];
 int    IJKPlusNeighbor[3];
 int    indexMinusNeighbor, indexPlusNeighbor;
 double valueMinusNeighbor, valuePlusNeighbor;

 // Loop over each coordinate 
 for (int i = 0; i< 3; i++)
  {
	 //Loop over all cell Ids
	for (int idC=0; idC<*(mesh_->NofCells()); idC++)
    {
		//compute the IJK index
		cellIDtoIJK(idC, &(IJKCurrCell[0]), cells);

		//compute the index of the current cell
		IJKMinusNeighbor[0] = IJKCurrCell[0];
		IJKMinusNeighbor[1] = IJKCurrCell[1];
	    IJKMinusNeighbor[2] = IJKCurrCell[2];

		IJKPlusNeighbor[0] = IJKCurrCell[0];
		IJKPlusNeighbor[1] = IJKCurrCell[1];
		IJKPlusNeighbor[2] = IJKCurrCell[2];

		// Adjust index to get the neighbours for the current coordinate
		IJKMinusNeighbor[i] = IJKCurrCell[i]-1;
		IJKPlusNeighbor[i] = IJKCurrCell[i]+1;

		// Apply periodic boundary conditions
		for (int c = 0; c < 3; c++)
		{
			if (IJKPlusNeighbor[c]>=cells[c])
			{
				IJKPlusNeighbor[c] = IJKPlusNeighbor[c] - cells[c];
			}
			if (IJKMinusNeighbor[c]<0)
			{
				IJKMinusNeighbor[c] = cells[c] + IJKMinusNeighbor[c];
			}
		}
    
		// Convert back to cell IDs
		indexMinusNeighbor = IJKtoCellID(IJKMinusNeighbor, cells);
		indexPlusNeighbor = IJKtoCellID(IJKPlusNeighbor, cells);

		//pull out values using the index
		valueMinusNeighbor = baseField_[indexMinusNeighbor];
		valuePlusNeighbor = baseField_[indexPlusNeighbor];

		//fill gradient fields
		if (i == 0)
		{
		*(gradx+idC) = (valuePlusNeighbor - valueMinusNeighbor)/(2*cell_size[i]);
		}
		else if (i == 1)
		{
		*(grady+idC) = (valuePlusNeighbor - valueMinusNeighbor)/(2*cell_size[i]);
		}
		else
		{
		*(gradz+idC) = (valuePlusNeighbor - valueMinusNeighbor)/(2*cell_size[i]);
		}
	}
  }

}
