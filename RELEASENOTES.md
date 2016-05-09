![logo](cpppo_logo.png)
======
A Compilation of Fluid-Particle Post Processing routines.

CPPPO is part of the [NanoSim Project](http://sintef.no/NanoSim)


Graz University of Technology and DCS Computing GmbH release CPPPO
------------------
31st of Decemmber 2014

The Graz University of Technology (TU Graz) together with DCS Computing GmbH (DCS) are please to announced the release of the `1.0.1-beta` version of the tool `CPPPO`. This version is distributed is licensed under the [Lesser General Public License](http://www.gnu.org/licenses/lgpl.html) by TU Graz and DCS.

###Features
Version 1.0.1-beta is the first public release of the CPPPO library and is meant to introduce its features to a wide audience of users and possible developers of this library. Specifically, the current version of CPPPO is able to

- compute a filtered (i.e., spatially-averaged) fluid velocity field during an OpenFOAM simulation run (i.e., “on the fly”) and in parallel. 
- read CSV data files (e.g., as dumped from ANSYS FLUENT) and perform filtering operations
- draw different types of samples (e.g., local fluid velocity, or two-point velocity correlations)
- each sample can be characterized with a number of markers

A number of test and verification cases, as well as documentation is supplied with CPPPO.

Note that CPPPO is still under active development, and the success of CPPPO relies to a large extend on contributions from users.

The developers of CPPPO are thankful for any comments, ideas, or contributions to currently available or future features of CPPPO that are posted in the CPPPO user forum](http://www.cfdem.com/forum).

###Important Notes
'beta' means, that this version of CPPPO has been tested internally to some extend, and now seeks for testing and application by external groups (i.e., non TUG or DCS). Thus, while TUG and DCS have carefully tested CPPPO, it is likely that users might find significant bugs, run into unexpected errors, or might miss information on how to use CPPPO. We are greatful for any users feedback, which should be directed to the [CPPPO user forum](http://www.cfdem.com/forum).

###Acknowledgement
Parts of the code were developed in the frame of the [NanoSim project](http://www.sintef.no/Projectweb/NanoSim/) funded by the European Commission through FP7 Grant agreement no. 604656.


