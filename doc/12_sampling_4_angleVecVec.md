Sampling command specifiers
======================


"type": "angleVecVec"
----------------------
"AngleVecVec" _sampling operations_ will compute the angle between two vectors for every cell, or between the slip velocity and the particle force for every particle. Sampled values are the module of the first vector, while the angle between the two vectors is used as _marker_. The user should consider that angles are given in degrees. Specific entries for "angleVecVec" are:
 
* "VecField1": requires a string of characters. The user defines here the fields to be sampled. At this point, the user can also use previously filtered fields. To access those fields, it is necessary to use the following syntax: (i) the name of the original field, and (ii) "\_", and finally (iii)  _command name_ corresponding to that _filtering operation_.

* "VecField2": requires a string of characters. The user can enter the name of an eulerian field (in this case all the eulerian domain will be sampled) or just enter "particelVelocity".With the latter option only a sample per particle will be drawn and the slip velocity module will be the sampled field. Filtered fields can be accessed (see "VecField1").

"AngleVecVec" can also be used to compute the angle between a vector field and the force acting on particle at particle centers. In this case, the user can choose if using the total force acting on particles or just a component. In this case, entries are:

 
* "VecField1": requires a string of characters. This will be the sampled field. At this point, the user can also use previously filtered fields. To access those fields, it is necessary to use the following syntax: (i) the name of the original field, and (ii) "\_", and finally (iii)  _command name_ corresponding to that _filtering operation_.

 
* "particleForce": requires a string of characters. If set to "TotalForce" the total force acting on particles will be used. This entry has to be set to "ForceModel" in order to sample with respect to a particular force model.

* "forceIndex": requires an integer. This field is necessary if "ForceModel" is used. The user has to enter the index that identifies the force to be used. Remember that the enumeration starts from zero.

Note that, if more than one vector field is declared in "VecField1", this _sampling operation_ will need as many _binning operations_  as the number of fields to be sampled (+1 if _formula_ is used) when "save2Bin" is 'true'. 

Example
-------
```
...
  "sampleA": {
                  "type"          : "angleVecVec",
                  "save2Disk"     : true,
                  "save2Bin"      : true,
                  "VecField1"     : "U_Favre U_Alg",
                  "VecField2"     : " U "
              
              },
              
   "sampleB": {
                  "type"          : "angleVecVec",
                  "save2Disk"     : true,
                  "save2Bin"      : true,
                  "VecField1"     : " U_Favre ",
                  "particleForce" : "ForceModel",
                  "forceIndex"    : 0
              
              },
...
```

Go back
-----------
 - [sampling](12_sampling.md)
