Selector command type
======================

The _selector command type_ defines selector objects in CPPPO. The corresponding _command specifiers_ can be _cell<type>_ or _filter_.
A few examples are:

* _cellIJK_ 

* _cellUnstruct_ 

* _filter_

* _cellRegion_

Cell selectors are detailed in their corrsponding documentation pages.

A command with the _filter specifier_ represents a _filter_. _filters_ have to be defined in the c3po.json file. For filters, cells are selected once their center is within the filter distance (specified in each direction, see below) from the master cell, i.e., the cell for which the filtered value is calculated. This means, that the filter width specified below is half of the box size over which the filtering (i.e., spatial averaging) is performed.

WARNING: It is recommended to define the filter width such that it extends a certain distance further than the cell center. In such a way, problems due to not selecting the cell can be avoided. Anyhow, CPPPO will add a small value (i.e., the value for `filterWidthTolerance` specified in the `mesh.json` file) to the filter width in order ensure the cell is selected in case the filter matches exactly the distance between cell centers.

Syntax
------
```
...

 "myFilter0":	  { 
                     "CoordSys": 0, 
                     "x":  1.0 , 
                     "y":  1.0 ,
                     "z":  1.0 ,
                     "selective":true,
                     "max":[10.5,  8.2,  5.3],
                     "min":[  -1,    4,  0.8]
                  },
...

```

Every _filter_ requires the following entry:

* "CoordSys": requires an integer value. If set to 0 the cell selection will be performed in cartesian coordinates and if set to 1, spherical coordinates will be used instead.

In the case of a cartesian coordinate system, three more entries are required:

*  "x": requires a double value. This entry specifies the spatial extension of the filter into the x direction.

*  "y": requires a double value. This entry specifies the spatial extension of the filter into the y direction.

*  "z": requires a double value. This entry specifies the spatial extension of the filter into the z direction.

In the case of a spherical coordinate system, one more entry is required:

*  "r": requires a double value. This entry specifies the spatial extension of the filter into the radial direction.

In addition, it is possible to just filter a box region out of the whole domain. In order to do that, three other entries must be added:
* "selective": requires a bool. If set to _true_ just a section of the domain will be filtered. 
* "max": requires an array of double values. The maximum box size.
* "min": requires an array of double values. The minimum box size.

For every _filter_ it is necessary to have a corresponding _cellIJK_ (or _cellUnstruct_) specifier. The corresponding _cellIJK_ (or _cellUnstruct_) specifier is selected based on the ordering in the `c3po.input` script. 



Example
-------

In the `c3po.input` file:

```
...

selector     cellUnstruct    cellSelector0
selector     cellIJK         cellSelector1

selector     filter          myFilter0
selector     filter          myFilter1 

...

```
While in the `c3po.json` file:

```
...

 "myFilter0":	  { 
                     "CoordSys": 0, 
                     "x":  1.0 , 
                     "y":  1.0 ,
                     "z":  1.0 
                  },
                  
 "myFilter1":	  { 
                     "CoordSys": 0, 
                     "x":  2.0 , 
                     "y":  1.0 ,
                     "z":  3.0,
                     "selective":true,
                     "max":[10.5,  8.2,  5.3],
                     "min":[  -1,    4,  0.8] 
                  },
...
```
In the previous example, _myFilter0_ will use _cellSelector0_ (the unstructured selector) and _myFilter1_ will use _cellSelector1_ (the IJK selector).
Note that the choice of the selector only depends from the mesh and the `c3po.input` should be:

```
...

selector     cellIJK         cellSelector0
selector     cellIJK         cellSelector1

selector     filter          myFilter0
selector     filter          myFilter1 

...

```         
In the case of a structured mesh without any refinement or a simple geometry. Alternatively it should have been:
```
...

selector     cellUnstruct         cellSelector0
selector     cellUnstruct         cellSelector1

selector     filter          myFilter0
selector     filter          myFilter1 

...

```      
In the case of an unstructured mesh or a complex geometry.

Go back
-----------
 - [commands](10_commandTypes.md) 
