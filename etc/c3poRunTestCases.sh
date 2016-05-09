#!/bin/bash
######################################################
# Shell script to automatically run some test cases
# of CPPPO
######################################################

source ./c3poFunctions.sh

printLogo


echo "Start CPPPO testing..."
echo ""

if  checkFile "../interface_CSV/c3po_csv" == "true" 
 then
 
 echo "Found stand-alone CPPPO-CSV interface"
 echo "will now run the test cases for the CSV interface"
 
 cd $C3PO_EXAMPLE_DIR
 cd ./codeTest
 cd csvTest
 ./Allrun.sh
 
 echo "did the case run correcty? - press enter to proceed."
 read
 
 cd $C3PO_EXAMPLE_DIR
 cd ./codeTest
 cd ./csvTestEulerData
 ./Allrun.sh
 
 echo "did the case run correcty? - press enter to proceed."
 read
 
else

 echo "ERROR: the CPPPO-CSV interface does not seem to be compiled"
 echo "please check or compile the library"
  

fi

cd $C3PO_EX_DIR

if checkEnv "$FOAM_USER_LIB" == "true" 

 then
 
 echo "Found OpenFOAM"
 echo "will now run the test cases for the OpenFOAM interface"
 
 cd $C3PO_EXAMPLE_DIR
 cd ./verification
 cd potentialStokesFilter
 ./Allrun.sh
 
 echo "did the case run correcty? - press enter to proceed."
 read
 
 cd $C3PO_EXAMPLE_DIR
 cd ./verification
 cd ./pisoFoam_c3po
 ./Allrun.sh
 
 echo "did the case run correcty? - press enter to proceed."
 read
 
else

 echo "ERROR: OpenFOAM does not seems to be installed"
 echo "cannot run test cases for the OpenFOAM interface"
  

fi

 cd $C3PO_SRC_DIR

