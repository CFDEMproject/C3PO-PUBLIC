#!/bin/bash

./Allclean
cp constant/transportProperties.stokes constant/transportProperties
blockMesh
snappyHexMesh -overwrite
decomposePar
mpirun -np 3 c3po_pisoFoam -parallel 2>&1 | tee log
mpirun -np 3 sample -parallel -latestTime

cd octave 
octave <  postproc.m
evince verificationAve.eps
evince verificationVar.eps
cd ..
