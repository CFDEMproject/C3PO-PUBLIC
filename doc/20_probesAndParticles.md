Registering probes and particles in CPPPO
======================

_CPPPO_ allows the user to define particles using the functions defined in the interface but it is also possible to define probes externally using a _json_ file. The main difference is that a probes is just defined by entering a position vector (i.e. they do not need the user to specify velocities, angular momenta or radii). In order to use probes or particles in _CPPPO_ it is necessary to set "useProbes" to _true_ in _c3po.json_ and edit the _probeSettings.json_ file in _c3po_control_. Additional _json_ files may be required as described below.

Syntax  
-------

In _c3po.json_

```
 "mainSettings":
   {
        "doFiltering": true,
        "doSampling":  false,
        "doBinning":   false,
        "dumpFormat": "hdf5",
        "verbose": false,
        "FieldsToRegister":
        {
        "Vectorfields": " U  ",
        "Scalarfields": "  "
        },
         "filterBCs":
        {
        "x":"periodic",
        "y":"periodic",
        "z":"periodic"
        },
        "useProbes" : true,
        "storageWriteFields"    :  false,
        "storageWriteParticles" :  true,
         "interfaceWriteFields" : false
     
   },
```
* "useProbes": requires a boolean value. Specifies if _CPPPO_ should look for _probeSettings.json_.

In _probeSettings.json_

```
{

 "probeSettings":
 {
 
  "probeNames": " samples particles",
  
  "samples":
  {
  
   "filters":"all",
   "readFromJson": true
  },
  
  "particles"
  {
   "filters":" filter0 filter1",
   "readFromJson": false
  
  }
 
 }

}


```
* "probeNames": requires a string of characters. The user has to specify the names used to identify different groups of probes/particles.
* The user has also to write an object for every particle group as shown in the example. Every object must contain the following entries:
 * "filters": requires a string. The user should enter the filters where the probe group is active.
 * "readFromJson": requires a bool. If set to _true_, _CPPPO_ will look for a file named _<groupName>.json_ were probe positions are detailed. If set to _false_ the user should register particles from the interface. 

Following the example, the user has to write another file for every probe group with "readFromJson" set to _true_.

In _samples.json_

```
{

  "probeData" :
  {
  "numberOfProbes":1,
  
  "positions":
  {
   "0": [ 4.001, 4.001, 4.001]
  }
 
 }



}

```
* "numberOfProbes": requires an integer. The user has to specify the number of probes to read.
* "positions": requires a list of samples. Samples need to be specified as follows:
 * "-sample_number-" : [ pos_x, pos_y, pos_z]
 
Samples have to be numbered starting from zero.

In C3PO_SRC/etc it is possible to find the _C3PO_createSamplingLocations_grid.m_ and _C3PO_createSamplingLocations_random.m_ to automatically generate the _filteringSamplingLocations.json_ file corresponding to an uniform distribution or a random distribution of sampling locations.


Example
-------

In _probeSettings.json_

```
{

 "probeSettings":
 {
 
  "probeNames": " samples particles",
  
  "samples":
  {
  
   "filters":"all",
   "readFromJson": true
  },
  
  "particles"
  {
   "filters":" filter0 filter1",
   "readFromJson": false
  
  }
 
 }

}
```

In _samples.json_

```
{

  "probeData" :
  {
  "numberOfProbes":1,
  
  "positions":
  {
   "0": [ 4.001, 4.001, 4.001]
  }
 
 }



}


```


Go back
-----------
 - [main](01_main.md) 

