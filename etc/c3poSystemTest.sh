#!/bin/bash

#===================================================================#
# sytsem settings test routine for c3po project 
#===================================================================#
source $C3PO_SRC_DIR/etc/c3poFunctions.sh


#- show gcc settings
checkGPP="true"

#- system settings
clear

printLogo

echo "**********************"
echo "CPPPO system settings:"
echo "**********************"

echo "CFDEM_VERSION=$CFDEM_VERSION"
echo "couple to OF_VERSION=$WM_PROJECT_VERSION"
echo "compile option=$WM_COMPILE_OPTION"

echo ""
echo "checking environment variables..."
checkEnvComment "$USEHDF5" '$USEHDF5' "yes"

echo ""
echo "checking if paths are set correctly..."
checkDirComment "$C3PO_SRC_DIR" '$C3PO_SRC_DIR' "yes"
checkDirComment "$C3PO_QT5_DIR" '$C3PO_QT5_DIR' "yes"
checkDirComment "$C3PO_HDF5_DIR" '$C3PO_HDF5_DIR' "no"
checkDirComment "$C3PO_HDF5_LIB" '$C3PO_HDF5_LIB' "no"
checkDirComment "$C3PO_CHEMKIN_SRC_DIR" '$C3PO_CHEMKIN_SRC_DIR' "no"
checkDirComment "$C3PO_EXAMPLE_DIR" '$C3PO_EXAMPLE_DIR' "yes"
checkDirComment "$C3PO_QT5_LIB" '$C3PO_QT5_LIB' "yes"
checkDirComment "$C3PO_QT5_INC" '$C3PO_QT5_INC' "yes"
checkDirComment "$C3PO_HDF5_INC" '$C3PO_HDF5_INC' "no"
checkDirComment "$C3PO_HDF5_INC" '$C3PO_HDF5_INC' "no"
echo ""

if checkDir "$C3PO_ADD_LIBS_DIR"
 then
 echo "using additional libraries defined in: "
 echo $C3PO_ADD_LIBS_DIR"/$C3PO_ADD_LIBS_NAME"
else
 echo "ERROR: it was not possible to find a valid definition of C3PO_ADD_LIBS_DIR"
 echo "please set this environment variable"
fi

echo "relevant library names:"
echo "*******************"

echo ""
echo "relevant aliases for c3po can be found here:"
echo $C3PO_SRC_DIR/etc/bashrc
echo "*******************"
echo ""

if [ $checkGPP == "true" ]
  then
    echo "g++:"
    which g++
    g++ --version

    echo "gcc:"
    which gcc
    gcc --version

    echo "mpic++:"
    which mpic++
    mpic++ --version

    echo "mpirun:"
    which mpirun
    mpirun --version
fi

