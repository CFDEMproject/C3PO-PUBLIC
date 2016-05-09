#!/bin/bash

#This script will run the c3po_pisoFoam application

./run.sh
./Allclean

echo " "
echo "-----------------------------------"
echo "average error at the last time step"
more error.txt
echo "-----------------------------------"

rm error*
