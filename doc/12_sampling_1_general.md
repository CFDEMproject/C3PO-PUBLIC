OPeration sampling general
======================


"type": "general"
----------------------

"General" _sampling operations_ will draw samples over all the Eulerian domain. The user can specify both the sampled field (which can be a vector or a scalar field) and the the _marker_ field (which has to be a scalar field).
The corresponding _command specifier_ is:

* _samplingGeneral_


Example
-------
In `c3po.input`:
```
...

operation samplingGeneral   sampleGeneral

...
```
In `c3po.json`:
```
...
                 
   "sampleGeneral": {
                  "marker"     : "interFace interFace_Favre",
                  "VFieldToSample": "U_Favre",
                  "SFieldToSample": " p ",
                  "component"   : 0,
                  "sampleCount" : -1,
                  "save2Disk"   : true,
                  "save2Bin"    : true,
                  "lagrangian"  :false
              },
...
```
Go back
-----------
 - [sampling](12_sampling.md)
