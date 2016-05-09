CPPPO-OpenFOAM interface
===============


Description
---------------------

The _OpenFOAM interface_ allows CPPPO functions to be called runtime during an OpenFOAM simulation. Furthermore, the _OpenFOAM interface_ is able to automatically read boundary conditions and to automatically register geometric fields (pass their addresses to CPPPO core). In order to couple CPPPO interface and OpenFOAM it is necessary to:

* Modify the ./c3po_control/mesh.json file

* Modify the OpenFOAM solver to include CPPPO functions

* Modifiy the test case (boudary conditions and domain decomposition have to be properly set)

It is important to keep in mind that, at the current state, CPPPO can only work with structured block-shaped meshes. Furthermore, all mesh elements have to be of identical size (i.e., no mesh refinement). and domain splitting among processors must produce subdomains of equal sizes. 

Modifications to mesh.json
--------------------------

The sub-fields below have to be added to the "mainSettings" field:

* "checkTolerance": requires a double value. Before running, CPPPO will check if all mesh volumes are equal. This field allows the user to set the tolerance.

Example
-------
```
{
"mainSettings": 
     {
        "checkTolerance":   1e-7,
        "verbose":          false
        
     }
}
```

Modifications to OpenFOAM solvers
---------------------------------

In order to CPPPO to run, the user has to create a _OpenFOAM interface_ object at the beginning of the solver. The user should also remember to destroy that object at the program end in order to completely free the heap memory. There are only two functions that the user needs to call:

* checkMyMesh(): this function has to be called once at the beginning before the solver loop.

* run(): this function starts a complete CPPPO run. It can be called multiple times.

If the user wants to register particle data in CPPPO, one more function has to be used.

* registerParticle(int particle\_id, double particle\_radius, double * particle\_position, double * particle\_velocity, double * particle\_force): This function will register one particle in CPPPO.

Furthermore, in order to link the _OpenFOAM interface_ to the solver the Make/options file in the solver folder has to be modified:

* EXE\_INC : the user has to add -I$(CPPPO\_SRC\_DIR)/core and  -I$(CPPPO\_SRC\_DIR)/interface\_OF to the include path.

* EXE\_LIBS : the user has to add -lCPPPO.

If HDF5 output is used, more additions are needed:

* EXE\_INC : the user has to add -I$(CPPPO\_HDF5\_DIR)/include.

* EXE\_LIBS : the user has to add -L$(CPPPO\_HDF5\_LIB), -lhdf5 and -lhdf5_cpp.

Go back
-----------
 - [main](01_main.md)
