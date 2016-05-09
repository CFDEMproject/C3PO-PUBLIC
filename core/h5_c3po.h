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

#ifndef C3PO_H5_NS_H
#define C3PO_H5_NS_H

#include <string>

namespace H5_C3PO_NS
{
   void createH5file(std::string filename_); 

   void createExtendibleDataset(std::string FILE_NAME,const char* datasetName_);

   void addOneArrayToH5(std::string FILE_NAME, const char* datasetName_, double* data, int NX);

   void addOneArrayToH5(std::string FILE_NAME, const char* datasetName_, int* data, int NX);

   void OneArrayToH5(std::string FILE_NAME, const char* datasetName_, double* data, int NX);

   void OneArrayToH5(std::string FILE_NAME, const char* datasetName_, int* data, int NX);

   void TwoArrayToH5(std::string FILE_NAME, const char* datasetName_, double** dataXY, int NX); 
   
   void TwoArrayToH5(std::string FILE_NAME, const char* datasetName_, double dataXY[][2], int NX); 

   void ThreeArrayToH5(std::string FILE_NAME, const char* datasetName_, double dataXY[][3], int NX); 

}
#endif
#endif
