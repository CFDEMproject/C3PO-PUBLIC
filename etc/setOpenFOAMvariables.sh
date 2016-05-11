#!/bin/sh
#----------------------------------*-sh-*--------------------------------------
# CPPPO
# Stefan Radl, Federico Municchi
# May 2016
#------------------------------------------------------------------------------
#
# Script
#     etc/setOpenFOAMvariables.sh
#------------------------------------------------------------------------------

#Additional libs required for Compilation of CPPPO with OF
if ! [ -z "$WM_PROJECT_VERSION" ] ; then
 if [ -z "$C3PO_ADD_LIBS_DIR" ] ; then
   if  [ -z "$CFDEM_ADD_LIBS_DIR" ] ; then
     export C3PO_ADD_LIBS_DIR=$C3PO_SRC_DIR/etc
   else
      export C3PO_ADD_LIBS_DIR=$CFDEM_ADD_LIBS_DIR
   fi
 else
    echo "using C3PO_ADD_LIBS_DIR=$C3PO_ADD_LIBS_DIR defined by user."
 fi



 if [ -z "$C3PO_ADD_LIBS_NAME" ] ; then
#if Additional libs does not exist, it has to be created
    if  [ -z "$CFDEM_ADD_LIBS_NAME" ] ; then

      if [ $WM_PROJECT_VERSION == "3.0.x" ]; then

       if [ -z "$USEHDF5" ]; then
        export C3PO_ADD_LIBS_NAME=additionalLibs_3.0.x
       else
        export C3PO_ADD_LIBS_NAME=additionalLibs_H5_3.0.x
       fi

     elif [ $WM_PROJECT_VERSION == "2.4.x" ]; then

        if [ -z "$USEHDF5" ]; then
         export C3PO_ADD_LIBS_NAME=additionalLibs_2.4.x
        else
         export C3PO_ADD_LIBS_NAME=additionalLibs_H5_2.4.x
        fi
     fi
    else
     export C3PO_ADD_LIBS_NAME=$CFDEM_ADD_LIBS_NAME  fi
    fi
 fi


# detect OF version


 if   [ -z $CFDEM_WM_PROJECT_VERSION  ] ; then

  if [ $WM_PROJECT_VERSION == "3.0.x" ]; then
    export CFDEM_WM_PROJECT_VERSION=30
  elif [ $WM_PROJECT_VERSION == "2.4.x" ]; then
    export CFDEM_WM_PROJECT_VERSION=24
  else
  echo "------------------------------------------------"
  echo "WARNING: Your OpenFOAM version is not valid"
  echo "         for CPPPO! Valid versions are:"
  echo "         3.0.x"
  echo "         2.4.x"
  echo "-------------------------------------------------"
  fi
 fi
fi
