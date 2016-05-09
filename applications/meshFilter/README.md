meshFilter
================

This is an application based on OpenFOAM to use the c3po library for filtering (i.e., spatial averaging) of data.

Insall
==========
- Make sure you have compiled CFDEM
- Go to `../CFDEMcoupling-RADL/src/c3po/`
- Type `./compileMe`
- To clean type `./cleanMe`


Example
========
- Example is located in `../CFDEMcoupling-RADL/src/c3po/examples/filterTest`
- Run `blockMesh` to build mesh
- You might want to change the mesh setup in `constant/polyMesh/blockMeshDict`
