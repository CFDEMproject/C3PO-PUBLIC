Radial Distribution Function Program Documentation “simpleRDF“
======================

Dated: 12/21/11
Goals
The program ‘simpleRDF’ calculates the radial distribution function based in a box (non-periodic, periodic, or with Lees-Edwards boundary conditions in one direction). The program code resides in ‘/home/shakti/sradl/masters/utilitiesOpenFOAM/simpleRDF’, and test cases can be found in ‘/home/shakti/sradl/masters/testCaseOpenFOAM/simpleRDF’. The code is structured as follows:
 
The folder ‘core’ contains relevant subroutines for the numerical operations of the program. In ‘dataTypes’ the files for the definition of data types and classes are saved. ‘IORoutines’ contains code for data exchange, and ‘Make’ is the directory that is needed by the make system (in our case we use ‘wmake’) to compile the code.
The main program resides in the file ‘simpleRDF.C’, the other *.h files are included in this file in order to keep the code tidy.
Main Algorithm
1.	Read a file/array of particle positions.
2.	Define the variables (box boundaries for all three dimensions, shear parameters)
3.	If the box is normal, use the RDF to find the distance between the particles.
4.	Execute rdfCore.h
Subroutines
Algorithm behind core/rdfCore.h
1.	Calculate the dimensions of the box in all 3 directions and store in the variable deltaBounds.
2.	Check for Lees-Edwards boundary conditions in all directions and calculate the corresponding deltaLEShift. 
3.	Subtract mirror images from deltaLEShift such that it is < deltaBounds: 
a.	For positive LEShifts, use the floor function. For negative LEShifts, use the ceiling function to correct for periodicity.
4.	Loop over the i, j particle pairs.
5.	Measure the distance between the particle pairs.
6.	Set the indicator for adding deltaLEShift to 0. Reset this indicator every time a different particle pair is being checked within the for loop.
7.	Correct for periodic boundary conditions in all directions.
a.	Add or subtract the box dimensions if the distance between the particle is smaller than or greater than the half the box dimensions.
b.	Set addDeltaLEShift to -1/+1 in each direction.
8.	Correct for LE boundary conditions.
9.	Check if the final distance is within the box dimensions.
10.	Bin the particle pairs. 
