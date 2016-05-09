# C3PO-PUBLIC
Public Version of C3PO Toolbox

![logo](cpppo_logo.png)
======
A Compilation of Fluid-Particle Post Processing routines.

CPPPO is part of the [NanoSim Project](http://sintef.no/NanoSim)


Copyright Notice
------------------

- Copyright 2014 - Graz University of Technology (S. Radl, F. Municchi).
- Copyright 2014 - DCS Computing GmbH, Linz (C. Goniva, C. Kloss).


All rights reserved.

License
-----------------
See the [LICENSE.md](LICENSE.md) file for details.

Warranty
-----------------
CPPPO is distributed in the hope that it will be applied and further developed by researchers and engineerings in academia or industry. However, CPPPO comes without any warranty, without even the implied warranty of merchantability or fitness for a particular purpose. 

Scope
---------------------------------------
CPPPO stands for “Compilation of Fluid-Particle Post Processing”, and is a publicly available library to analyze simulation data on the fly, e.g., spatially average (i.e., "filter") fluid velocity fields, draw samples of forces on particle, and compute running averages (i.e., "bin the data"). For example, CPPPO can compute filtered velocities fields and can computed advanced statistical information from this filtered data (e.g., filtered drag coefficients).

CPPPO is available via https://github.com/CFDEMproject, and is designed as an add-on to OpenFOAM and/or CFDEM-based simulations. Specifically, CPPPO is designed as a 'core' library, with various interface libraries that are meant to be instantiated in your simulation software. As such, CPPPO is part of the [CFDEMproject](http://www.cfdem.com) initiated by Christoph Kloss and Christoph Goniva. However, CPPPO can also be run in stand-alone mode for data analysis, e.g., from ANSYS FLUENT output.

Hints for Usage
-----------------

- the 'doc' folder contains information on how to use CPPPO, and the 'examples' folder contains test cases that illustrate the usage of CPPPO.
- for installation instructions see INSTALL.md
- there are separate README.md files in each sub-directory, so you might want to check them
- applications are available in applications/
- in case you plan to add more details to the documentation, Markdown should be used for this purpose. The tool "ReText" (in the OpenSUSE standard package) or gedit (with the appropriate [gedit-markdown plugin](http://www.jpfleury.net/en/software/gedit-markdown.php)) for previewing markdown documentation on linux systems should be used to edit and view documentation.

Credits
-------------------
The following persons have made significant contributions to CPPPO (in alphabetic order)

- Christoph Goniva (DCS)
- Federico Municchi (TU Graz)
- Stefan Radl (TU Graz)
