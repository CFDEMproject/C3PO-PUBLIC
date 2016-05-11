 Installation Instructions
==================

Linux Environment
-----------------------
Be sure to have all the required libraries installed before compiling CPPPO.
You will need, at least, Qt-5 and an MPI library in order to compile CPPPO

### Qt (required)
Be sure you have the Qt library installed:

> Visit [QT Download Page](http://www.sysads.co.uk/2014/05/install-qt-5-3-ubuntu-14-04/)

> Type (e.g. for 64 bit opperating system) `wget http://download.qt-project.org/official_releases/qt/5.3/5.3.0/qt-opensource-linux-x64-5.3.0.run`

> Type `chmod +x qt-opensource-linux-x64-5.3.0.run`

> Type `./qt-opensource-linux-x64-5.3.0.run`

### openMPI (required)
Also Intel MPI libraries can be used in place of openMPI.
OpenMPI is surely available in your system repositories, we
will just present the installation for Ubuntu and OpenSUSE.

Ubuntu: type `sudo apt-get install libopenmpi-dev openmpi-bin`

OpenSUSE: type `sudo zypper install openmpi-devel`

CPPPO was shown to work with openmpi-1.8.8.

### HDF5 (optional)
In case you want to use HDF5 capability, be sure it is installed:

> Check if the HDF5 libraries and development files are available in your package manager. Be sure you have version 1.8.15 or higher available (note: earlier versions might not work!). Install them as root (preferred option).

In case your package manager does not provide and appropriate version of HDF5, please visit [HDF5 Download Page](http://www.hdfgroup.org/ftp/HDF5/current/src/unpacked/release_docs/INSTALL)

> Download the last released version in your $HOME directory

> Unpack typing `gunzip < hdf5-X.Y.Z.tar.gz | tar xf -` where X,Y and Z are the package numbers

> Type `cd hdf5-X.Y.Z`

> Configure the library using `./configure --prefix=$HOME/hdf5-X.Y.Z --enable-cxx`

> Type `make`

> Type `make install`

WARNING: Be sure you have HDF5 version 1.8.15 or higher installed!

Most important, in your .bashrc file, please include a section to set the main switch for including HDF5 libraries:

> unset  USEHDF5

> export USEHDF5=true

Please make sure that an "additionalLibs" file is found and adapted to your system by defining the environmant variable CFDEM_ADD_LIBS_DIR. It is used to link additional libraries to the applications.

Also, it is useful to have a HDF5 viewer installed. Visit [this homepage](http://www.hdfgroup.org/products/java/release/download.html) to install 'hdfview'. Just download the installation script, and follow the installation instructions.

### Octave & JSONLAB (optional)
We recommend using octave version 3.8.x or later for correct display and printing of result graphs.

Be sure you have correctly set up Octave (including `JSONLAB`) for post processing

> Visit [jsonlab Download Page](http://sourceforge.net/projects/iso2mesh/files/jsonlab/)

> Download and unpack latest Version (here we assume that you save it to `'/home/username/utilities/jsonlab'`)

> Go to your $HOME folder and type `touch .octaverc`

> Edit `.octaverc` (e.g. `vim .octaverc`) and add the path to jsonlab folder by typing e.g. `addpath('/home/username/utilities/jsonlab')`

> Save `.octaverc`, close and countercheck the added path by typing `octave`

> In Octave type `path` to check all loaded paths - `/home/username/utilities/jsonlab` should appear on top

### Environment Variables
Be sure you have correctly set the variables C3PO_SRC_DIR, e.g., in your .bashrc you should have:

>export C3PO_SRC_DIR=$HOME/C3PO-PUBLIC  

>export C3PO_ADD_LIBS_DIR=$C3PO_SRC_DIR/etc  

>export C3PO_ADD_LIBS_NAME=additionalLibs_3.0.x 

>export C3PO_QT5_DIR=$HOME/Qt/5.3/gcc_64 

>export C3PO_QT5_LIB=$C3PO_QT5_DIR/lib 

>export C3PO_QT5_INC=$C3PO_QT5_DIR/include 

>export C3PO_HDF5_DIR='pathToTheHDF5InstallationDir' 

>export C3PO_HDF5_LIB=$C3PO_HDF5_DIR/lib 

>export C3PO_HDF5_INC=$C3PO_HDF5_DIR/include 

>export USEHFD5=true 

>export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$C3PO_QT5_LIB:$C3PO_HDF5_LIB 

>. $C3PO_SRC_DIR/etc/bashrc

Remember to check the environment variables carefully since they may be different in your system!

In case OpenFOAM(R) is not linked, the following variables are not necessary:
>C3PO_ADD_LIBS_DIR

>C3PO_ADD_LIBS_NAME


In case the HDF5 library is not linked, the following variables are not necessary:
>C3PO_HDF5_DIR

>C3PO_HDF5_LIB

>C3PO_HDF5_INC

>USEHFD5


The C3PO_ADD_LIBS_DIR must point to an existing directory of the user's choice. The directory must contain a file called "additionalLibs" that is used to defined which libraries are linked the OpenFOAM(R)-type applications in the CPPPO package. The user has to provide this file
when OpenFOAM(R) applications are linked.
The user should also set C3PO_ADD_LIBS_NAME accordingly to its additionalLibs file.
Template are available in etc/ .

WARNING: In case the compiler can not find the HDF5 library path try replacing /lib with /lib64 in C3PO_HDF5_LIB:

>export C3PO_HDF5_LIB=$C3PO_HDF5_DIR/lib64

Note that we just described standard paths to /include and /lib folders. Paths in your system may be different.

The lines 'export C3PO_HDF5_DIR=' is used to decide whether CPPPO is compiled with or without HDF5 capability: if C3PO_HDF5_DIR is not set, C3PO will not compile its HDF5 modules, and hence will not link a HDF5 library (see the 'compileMe' script!). To force that the HDF5 modules are NOT compiled, you can set the environmental variable 'USEHDF5' to false (see the HDF section above).

In case your compiler cannot find mpi.h, set the variable

>MPI_INCLUDE_PATH

in your .bashrc such that it points to the include patch of your OpenMPI installation.

General Hints
------------------------------------
the C3PO packages consists of a core library (in ./core), interface modules that are compiled as a library (e.g., ./interface_OF), as well as sample applications (in ./applications).

To build all libraries and applications, use the

> c3poComp

script. Note, that the 'compileMe' script also triggers the installation or de-installation of the HDF5 modules. This script will also generate the Makefile.lib file.

> c3poClean

alias will call the cleanMe script and clean all libraries and applications, as well as reset the make files to NOT include the HDF5 modules.


Building the C3PO Core Main Application ./core
------------------------------------
The C3PO core modules will be build as an executable with just one instance of the C3PO core library.
Currently, only one linux architecture is supported.  To compile, simply do:

>make fedora_fpic

Note, a

>make clean-all

is recommended to clean out pre-compiled files.

Building the C3PO Core Library in ./core
------------------------------------
C3PO will be build as library. Currently, only one linux architecture is supported. To compile, simply do:

>make makelib

>make - Makefile.lib fedora_fpic

Note, a

>make clean-all

is recommended to clean out ALL (including third party libraries) pre-compiled files. Note, that the C3PO core library will be build as a static library, and then linked to the interface library (e.g., the interface_OF library).

You might want to set a symboli link to the main application of c3po core using "ln -s 'TARGET' 'LINK_NAME'".


Building the C3PO OpenFOAM Interface libary in ./interface_OF
------------------------------------
The interface library will be build using the OpenFOAM build system

>wmake libso

and cleaned with

>wclean

Be sure you freshly compile  the OpenFOAM library, i.e., always do a "wclean" before "wmake libso"!

Building the C3PO Test Applications in ./applications
------------------------------------
OpenFOAM-style applications can be build using

>wmake

in the respective directory
