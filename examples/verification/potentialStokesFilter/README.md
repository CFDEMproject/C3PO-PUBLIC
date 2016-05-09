POTENTIAL-STOKES verification case for CPPPO
===============


Description
---------------------
This test case will impose the Stokes flow around a shpere on a cubic domain and CPPPO will compute the Favre averaged velocity field and its Favre Variance at the particle center.
Results are plotted using octave against analytical results. The procedure is repeated when an irrotational flow is imposed instead.

Usage
---------------------
The verification case can be run using the _Allrun.sh_ script (which will run both the flow fields automatically cleaning the case at the end) or using the _runStokes.sh_ and _runPotential.sh_ scripts to run just one of the cases (the case will not be automatically cleaned at the end).
