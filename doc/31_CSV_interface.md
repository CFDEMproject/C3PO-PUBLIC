CPPPO-CSV interface
===============


Description
---------------------
The _CPPPO CSV interface_ allows the user to "link" CPPPO core to every simulator capable of csv data output. When compiled, this interface module will be packaged into a stand-alone application (see the source code in CPPPO/interface\_CSV/main.cpp). This application will automatically read and process csv data files from a folder named "data". Those csv files MUST contain cell center coordinates and every vector field has to printed with adjacent components. The user has also to provide CPPPO with adequate information on mesh and cell size. The user should also note that _CPPPO CSV interface_ can run in parallel. In this case CPPPO will divide the domain in slices.

The _CPPPO CSV interface_ requires two or three CSV files:

* "data" folder : CPPPO will look for a folder named "data" and will read all the .csv files inside it in alphabetic order. CPPPO will generate a list named "fileList.c3po" that will be used to read the data sequentially. Every file requires an header with the field names. Vector field components have to be in adjacent columns, and the x-vector must contain the string '-x'.
* "mesh.csv" : this file contains mesh data and it requires at least four columns with an header (_CPPPO_ will always skip the first row). Three columns have to contain x,y and z coordinates of cell centres while one column has to contain the cell volumes. The user can manually specify the columns to read editing the "CsvColumns" field in "mesh.json"; if this field is not present, CPPPO wil automatically read the first four columns with the following order: x, y ,z ,volume.
* "lagrangian.csv" : this file is optional and contains particle data. It has to be structured as follows:
 * the first column holds the particle radii.
 * columns from two to four hold the particle position in coordinates.
 * columns from five to seven hold the particle velocity components.
 * columns from eight to ten hold the components of the torque acting on particles.
 * the remaining columns can hold an arbitrary number of forces acting on particles. For every force, the three components have to be provided.

NOTE: cell indexes are not necessary and the input file should just contain field data. The _CPPPO CSV interface_  can process data from two-dimensional simulations but the mesh.csv file has to contain the three coordinates and
the volumes anyway. Field data can be provided only for two components but the output will contain an additional component with all the calues set to zero. No one-dimensional data are allowed. 
Syntax
---------------------
The _CPPPO CSV interface_ require mesh parameters to be specified in mesh.json "mainSettings":

* "MaxDomain": requires an array of doubles. The user has to enter the x, y and z domain max coordinates.

* "MinDomain": requires an array of doubles. The user has to enter the x, y and z domain min coordinates.

* "domain decomposition": requires an integer. The user can enter '0', '1', '2' corresponding to x, y and z direction. If run in parallel, CPPPO will decompose the domain over the specified direction. Deault is '2'

* "Cell Size": requires an array of positive doubles. The user has to enter the x, y and z cell dimensions. 

* "CsvColumns": optional, requires four sub-entries. This field can be edited to specify the columns CPPPO has to read. Consider that the first column is 0. Sub-entries are:
 * "cellx": requires an integer. The column corresponding to the x-coordinate of the cell centres. 
 * "celly": requires an integer. The column corresponding to the y-coordinate of the cell centres. 
 * "cellz": requires an integer. The column corresponding to the z-coordinate of the cell centres. 
 * "cellV": requires an integer. The column corresponding to the cell volumes. 



Example
---------------------

From `c3po.json`:

```
"mainSettings":
   {
        "doFiltering": true,
        "doSampling":  true,
        "doBinning":   true,
        "dumpFormat": "hdf5",
       "FieldsToRegister":
        {
        "Vectorfields": "   x-coordinate ",
        "Scalarfields": " pressure solids-vof"
        },
        "verbose": true
   },
   
...
```

From `mesh.json`:

```
{
   "mainSettings": 
     {
        "checkTolerance":   1e-7,
        "verbose":          false,
        "MaxDomain": [ 0.12, 0.12, 0.12],
        "MinDomain": [ -0.12, -0.12, -0.12],
        "Cell Size": [ 0.01, 0.01, 0.01],
        "domain decomposition": 0
        "CsvColumns":
        {
         "cellx" : 2,
         "celly" : 1,
         "cellz" : 3,
         "cellV" : 4
        }
     }


}
```

Go back
-----------
 - [main](01_main.md)
