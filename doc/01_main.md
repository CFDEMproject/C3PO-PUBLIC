CPPPO Documentation
=====================

The CPPPO (or in short "c3po") package is build on the philosophy of on-line (on the fly) filtering of data. CPPPO is packaged in a form such that i can be easily linked with different simulators for fluid and particle flow. The key novelty of CPPPO is that all data operations (e.g., data filtering) can be specified at run time, i.e., no expert knowledge (from a coding point of view) is necessary to gather and process averaged or filtered data. Furthermore, input data can be provided by means of .csv files allowing CPPPO to interface with virtually every data that is generated in various simulators.
Currently supported simulators are:

* OpenFOAM (via a dedicated interface)
* ANSYS FLUENT (via data files in .csv format)

The CPPPO source code is split into three major parts: (i) the core, (ii) interface modules, and (iii) example applications.

CPPPO core
-------------

This is the core library that can be linked to an interface module of choice. The core library provides functionality for 

* Selecting cells or particles to be averaged over (so-called "selectors").
* Spatially-anisotropic filtering on a structured/unsitructured grid.
* Calculation of gradients and shear rates of filtered fields.
* Operations on solid and fluid domain  (i.e., "operations", such as (i) filtering, (ii) sampling, and (iii) binning of data).
* Input and output (to screen, log files, as well as JSON and HDF5 files).
* MPI data exchange. CPPPO mainly uses collective communications in order to make data available to all processors during data filtering operations.
 
CPPPO core also contains a standalone application, which can be used to perform unit tests.

CPPPO Interface Modules
----------------

These libraries contain references to the simulator's data, and are the driver to the functionality provided by CPPPO core. The interface modules feed data into CPPPO core, and fill received data into appropriate containers that can be saved by the simulator. Two interfaces are currently available:

* OpenFOAM interface (`interface_OF`): provides functions to read boundary conditions, mesh and cell size, passing or creating geometric fields to CPPPO and run CPPPO from users applications.
* CSV interface (`interface_CSV`): provides functions to read and import to CPPPO data in .csv format as well as .csv output. A parser function allows to order .csv data in CPPPO format. Mesh information has to be provided by the user, though. This interface module comes as a standalone application.

CPPPO Applications
-------------------

CPPPO provides applications that highlight how CPPPO Interface Modules can be integrated into simulators (they all require OpenFOAM):

* MeshFilter: This application provides an example of CPPPO usage for filtering eulerian fields.
* CPPPOpisoFoamScalar: This application shows how to link CPPPO to a full OpenFOAM solver.
* pisoFoam_c3po : Similar to CPPPOpisoFoamScalar but without the scalar transport.
* postProCFDEM : Application to post process CFDEM data (requires CFDEMCoupling).
* irrotationalFilter : Imposes a potential flow around a sphere and launch CPPPO.
* stokesFilter : Imposes a Stokes flow around a sphere and launch CPPPO

Main index
----------

|  Basics                              | Operations                    | Selectors  | Interfaces                                      | Advanced  |
|:---                                  |:---                           |:---        |:---                                             |:---       |
| [Installation](INSTALL.md)           | [Filtering](11_filtering.md)  | [Selector type](14_selector_type.md)            | [OpenFOAM interface](30_OpenFOAM_interface.md)  |[Gradient calculation](03_postFilteringOperationsGradients.md) |
|  [Input files](02_c3poInput.md)      | [Sampling](12_sampling.md)   | [Cell selector](14_cellSelector.md)           | [CSV interface](31_CSV_interface.md)            |[Probes and particles](20_probesAndParticles.md) |
|  [Command types](10_commandTypes.md) | [Binning](12_binning.md)    | [Region selector](14_regionSelector.md)             |                                                 |                                                 |

 
Go back
-----------
 - [README](00_README.md)
 
