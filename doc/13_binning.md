Binning command specifiers
======================


_Binning command specifiers_always need to be followed by the _command name_. Any _operation_ declared using a _binning command specifier_ is called _binning operation_ and has to be defined in the `c3po.json` file. _Binning operations_ perform running statistics using data provided by _sampling operations_ after those are divided in groups called _bins_. _Binning operations_ calculate:

* the number of samples in each bin

* the running mean value

* the running variance

Syntax  
-------
```
...

 "bin2PC": {
                  "bincount"        : 180,
                  "binUpBorder"     : 10.2,
                  "binLowBorder"    : -1.3,
                  "overwrite"       : false
                  
             },
 
...
```

The number of entries needed to define a _binning operation_ are:

* "bincount":       requires an integer value. This is the number of intervals in which sampled values are divided. 

* "binUpBorder":    requires a double value. This value represents the maximum value of the _marker_ accepted by this _binning operation_.

* "binLowBorder":   requires a double value . This value represents the minimum value of the _marker_ accepted by this _binning operation_. 

* "overwrite": requires a boolean value. If set to true, all the files dumped by this _binning operation_ will be overwritten the next time the bins are supplied with data (e.g., when a sampling operation is called). Furthermore, old bins will not be cleared, resulting in a statistics that are time-averaged over the whole simulation time. If it is set to false, all the files dumped by this _binning operation_ will NOT be overwritten. Also, the statistics recorded by the _binning operation_ will be cleared, resulting in multiple data files with statistics that are representative for on time step only. 
If more than one marker is used (i.e this binning operation is connected to a multimarking sampling operation) the following syntax has to be adopted:

```
...
"bin": {
              "NumOfMarkers":2,
              
                "Marker0" :
                 {    
                  "bincount": 20,
                  "binUpBorder": 1.5,
                  "binLowBorder": 0
                  
                  },
               
                
                "Marker1" : 
                { 
                  "bincount": 5,
                  "binUpBorder": 1.1,
                  "binLowBorder": 0.8
                  
                },
                  
                  
              "overwrite": false    
             },
...
```

The additional entries are:

* "NumOfMarkers":       requires an integer value. The number of markers used for this _binning operation_. 

* "MarkerX":            where X is the marker id (integer, starts from zero). The marker id is not arbitrary; markers have to be ordered (as in the example).
                        Every "MarkerX" requires a set of sub fields corresponding to those used in the case of a single marker (as shown in the example).
                        
Output  
-------
Output from _binning operations_ can be found in the _c3po\_binning_ folder. Center values and global running statistics are dumped to disk. In the case "overwrite": false a number indicating the time step is appended at the end of file. When multimarking is used, every file will contain the global statistics with respect to the first marker for a particular combination of values from other markers. The user has to imagine that theese files refer to a parameters space and are indexed in an IJK fashion (even if, of course, the parameters space can have more than three dimensions). Using this notation, the output will be in the form of :

```
bin_proc{ processor number }_{ filterName }_id_{ binOfMarkerN binOfMarker(N-1) ... binOfMarker1 }_global{ numberOfTimeStep (if "overwrite": false) }.h5
```
To avoid CPPPO to dump an excessive number of files, the user should define "Marker0" the one with the highest bincount.

Go back
-----------
 - [commands](10_commandTypes.md) 
