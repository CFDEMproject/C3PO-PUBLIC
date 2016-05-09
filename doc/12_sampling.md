Sampling command specifiers
======================


_Sampling command specifiers_ always need to be followed by the _command name_ . Any _operation_ declared using a _sampling command specifier_ is called _sampling operation_ and has to be defined in the `c3po.json` file.
There are several _Sampling command specifiers_ according to the operation type.
However all the sampling operations share some similar entries in `c3po.json`.

Syntax  
-------
```
...

 "sample2PC": {
                  "marker1"     : "interFace",
                  "VFieldToSample": "U_Favre",
                  "SFieldToSample": "Y",
                  "component" : 0,
                  "samplesLimiter"  : [0,20],
                  "sampleCount" : 20,
                  "save2Disk"   : true,
                  "save2Bin"    : true,
                  "lagrangian"  : false
              },
...
```
Every _sampling operation_ requires the user to specify some fields in the `c3po.json` file according to the "type" entry. However, any _sampling operation_ requires the following fields (for more details on type-specific entries see the corresponding *.md file):

* "save2Disk": requires a boolean value. If set to true, the sampled values will be written to disk in hdf5 or json format according to the specification provided in the "mainSettings" section of the `c3po.json` file. The corresponding output will consist of two arrays: the first is a marker, and the second one is the corresponding sampled value.

* "save2Bin": requires a boolean value. If set to true, the sampled values will be sent to the corresponding _binning operation_. If you choose this option, please be sure that a corresponding _binning operation_ exists. Otherwise, an error will be thrown.

* "VFieldToSample": requires a string of characters. The user has to specify the vector fields he/she wants to sample. At this point, the user can also use previously filtered fields. To access those fields, it is necessary to use the following syntax: (i) the name of the original field, and (ii) "\_", and finally (iii)  _command name_ corresponding to that _filtering operation_. 

* "component": requires an integer value. This field is only required if one or more vector fields are sampled. If set to 0, the x component will be sampled. If set to 1, the y component will be sampled. If set to 2, the z component will be sampled.

* "SFieldToSample": requires a string of characters. The user has to specify the scalar fields he/she wants to sample. At this point, the user can also use previously filtered fields. To access those fields, it is necessary to use the following syntax: (i) the name of the original field, and (ii) "\_", and finally (iii)  _command name_ corresponding to that _filtering operation_. 

* "marker": requires a string of characters. The user can define which field data is used to mark each sample (to be used for calculating conditional averages). Markers can also be previously filtered data fields, which can be selected similarly to the "fieldToSample" specification.

* "sampleCount": requires an integer. The user can specify the number of cells to sample. If set to -1, the whole domain will be sampled.

* "samplesLimiter": requires an array of two doubles. The user enter simple limiters (first value is the minimum and second value is the maximum) for the sampled values (for more complecated limiters, refer to the _formula_ documentation).

Note that each _sampling operation_ requires as many _binning operations_ as the number of fields to be sampled ( +1 if a _formula_ is used) when "samve2Bin" is set to 'true'.
This sampling operation allows multimarking in order to calculate conditional expected values of the sampled fields in an arbitrary parameter space. Results are dumped to disk in "c3po_binning" as explained in "13_binning.md". In addition a multimarking _sampling operation_ must be connected to the corresponding multimarking _binning operation_.  

It is also possible to sample inside a specific region of the domain adding the following entries:
* "selective": requires a bool. If set to _true_ just a section of the domain will be sampled. 
* "max": requires an array of double values. The maximum box size.
* "min": requires an array of double values. The minimum box size.
 

Note, that documentation on how to use the different sampling operations delivered with CPPPO is provided in separate files.
Some sampling operations will allow the following entry:
* "lagrangian":           requires a boolean value. If set to _false_, every cell in the domain will be sampled, otherwise only cells at particle/probe centres will be sampled. If this is set to _true_ the name of the probes/particles group should be provided with an additional entry:
 * "probesName": requires a string of characters. The user has to enter the probes/particles group that will be used by this operation. Note that just one group per operation can be defined.


In addition, each sampling operation allows to specify a formula parser as explained in  "[using the formula parser](12_sampling_5_formula.md)".

Example
-------
```
...
                 
   "sample2PC": {
                  "marker"     : "interFace",
                  "vectorFieldToSample": "U_Favre",
                  "component" : 0,
                  "sampleCount" : 20,
                  "sampleDelta" : 0.5,
                  "save2Disk"   : true,
                  "save2Bin"    : true,
                  "lagrangian"  : false
               },
  
 
                 },
...
```

List of _Sampling command specifiers_
-----------
 - [general sampling](12_sampling_1_general.md)
 - [particle sampling](12_sampling_2_particle.md)
 - [two-points correlation sampling](12_sampling_3_twoPointCorr.md)
 - [angular sampling](12_sampling_4_angleVecVec.md)
 
 
Go back
-----------
 - [commands](10_commandTypes.md) 
