#!/bin/bash
mpirun -np 3 c3po_pisoFoam -parallel 2>&1 | tee log
