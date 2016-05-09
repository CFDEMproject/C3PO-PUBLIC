#!/bin/bash

#This script will run the stokesFilter application
./Allclean

blockMesh
setFields
#run in parallel only when all selectors are cellUnstruct. This because the decomposition produces an unstructured mesh.
#decomposePar -force
#mpirun -np 3 irrotationalFilter -parallel 
#snappyHexMesh -overwrite
stokesFilter
cd octave
cp postproc.m.stokes postproc.m
#substiture 'processor0' with 'processor1' in postproc.m when run in parallel
octave postproc.m

