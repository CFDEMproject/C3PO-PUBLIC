Region Selectors
======================

 A _cellRegion_ selector enables the user to analyze irregular region within the computational domain. The mesh must use IJK organized cells. The region will be typically selected based on an unfiltered scalar field (e.g., the voidfraction field) and a threshold value. The only drawback of this selector is the speed, since a thresholding is required. The user can decide, however, if the regions are identified every time step, or just once (i.e., based on the first time step). Also, currently this feature is limited to single-core runs.

* _cellRegion_

Cell region selectors only work for cartesian grids!
Cell region selectors do not require any _filter_ to run.


Syntax
------

```
...
 
 "myBubbles":
   {
    "fieldForSelection":"voidfraction",
    "threshold":0.5,
    "identifyBelowThreshold": false,
    "updateRegionsEachStep": false
   
   }

...

```

Each region selector requires some entries in `c3po.json`:

* "fieldForSelection": requires a string of characters. It is the field used to check if the cell is inside the region.
* "threshold": requires a double. The numerical value of the threshold.
* "identifyBelowThreshold": requires a bool. If set to true it will consider inside the region each cell with a value of the "fieldForSelection" below "threshold".
*  "updateRegionsEachStep": requires a bool. If set to false it will check for the region just the first run.

Output
------

Output from _cellRegion_ selectors is available in the `c3po_dataStorage_regions` folder.
Each region produces two output files:

* \[regionName]\_\[procNumber]\_cellList : contains a list of cells for the corresponding processor.
* \[regionName]\_\[procNumber] : contains region data like the center of mass, the total volume and the toatal surface.


Example
-------

In the `c3po.input` file:

```
...

selector  bubble        myBubbles

...

```
While in the `c3po.json` file:

```
...
"myBubbles":
   {
    "fieldForSelection":"voidfraction",
    "threshold":0.5,
    "identifyBelowThreshold": false,
    "updateRegionsEachStep": false
   
   }

...
```


Go back
-----------
 - [main](01_main.md) 
