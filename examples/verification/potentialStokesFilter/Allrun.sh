#!/bin/bash

#This script will run the stokesFilter application
./Allclean
./runStokes.sh
./Allclean
./runPotential.sh
./Allclean

echo " "
echo "-----------------------------"
echo "average error for Stokes flow"
more errorStokes.txt
echo "-----------------------------"
echo "average error for irrotational flow"
more errorPotential.txt
echo "-----------------------------"

rm error*
