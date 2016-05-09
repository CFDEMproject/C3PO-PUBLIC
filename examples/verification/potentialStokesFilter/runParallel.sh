#!/bin/bash

#This script will run the stokesFilter application
./Allclean

blockMesh
setFields
decomposePar 
#run in parallel only when all selectors are cellUnstruct. This because the decomposition produces an unstructured mesh.
#decomposePar -force
#mpirun -np 3 irrotationalFilter -parallel 
#snappyHexMesh -overwrite
mpirun -np 3 stokesFilter -parallel
cd octave
cp postproc.m.stokes postproc.m
#substiture 'processor0' with 'processor1' in postproc.m when run in parallel
octave postproc.m
evince verificationAve.eps
evince verificationVar.eps

