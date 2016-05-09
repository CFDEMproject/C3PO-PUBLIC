Cell Selectors
======================

Cell selectors allow to select a collection of cells in the computational domain according to the mesh topology.
They corresponding _command specifiers_ are: 

* _cellIJK_ 

* _cellUnstruct_ 



A command with the _cellIJK_ specifier represents a cell selector based on IJK indexing. A command with the _cellUnstruct_ specifier represents a cell selector based on unstructured mesh. There are significant differences and
usages between the twos:

* A _cellIJK_ selector is based on the assumption that all the cells are equal (same volumes, same lengths) and hexahedrals. Furthermore, the domain is hexahedral itself. If this is the case, any cell can be identifyed using an IJK notation. This results in a fast algorithms because, if a particular region in the domain is given, the number and ids of the cells inside that region are calculated using simple relations.

* A _cellUnstruct_ selector does not pose any limitation on the shape and lenghts of cells and domain. The only drawback of this selector is the speed which is lower compared to the _cellIJK_ selector when the filtersize is small. 

Both these selectors require the domain to be properly decomposed. Sub-domains can have any shape (if the _cellUnstruct_ seector is used) but the hexahedrals built using maxima and minima of every subdomain should not overlap. In other words, processor-processor boundaries should be perpendicular or parallel to the main axes. 

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
 - [main](01_main.md) 
