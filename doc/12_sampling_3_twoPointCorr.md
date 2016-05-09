Sampling command specifiers
======================

"type": "twoPointsCorr"
----------------------
"twoPointsCorr" _sampling operation_s will compute numerators and denominators for calculating the trace of a two-point correlation tensor of a vector field. Additional entries for "twoPointsCorr" _sampling operations_ are:

*  "fieldToSample": requires a string of characters. The user has to specify the VECTOR field he wants to sample. At this point, the user can also use previously filtered fields. To access those fields, it is necessary to use the following syntax: (i) the name of the original field, and (ii) "\_", and finally (iii)  _command name_ corresponding to that _filtering operation_.

*  "filteredField": requires a string of characters. If the user wants CPPPO to compute the "twoPointsCorr" of the fluctuating component of a vector field (i.e. the field value minus the filtered field value), he/she should fill this field with the name of the corresponding filtered field. To call a filtered field, it is necessary to use the syntax described in the "fieldToSample" specification.

* "phaseFractionField": requires a string of characters. The user has to specify, in this field, the name of the the scalar field representing the phase fraction to identify if a certain sample is valid (e.g., if calculating a 2-point correlation of a particle quantity, the sample can be only drawn at locations where particles are present). The "phaseFractionField" can also be a previously filtered data field, which can be selected similarly to the "fieldToSample" specification.

* "sampleCount": requires an integer value. This field represents the number of samples to take in every direction.

* "sampleDelta":requires a double value. The user has to specify, in this field, the distance (in terms of the same unit lenght used in the mesh) between two samples.

* "lagrangian": requires a boolean value. If set to true the sampling will start from every cell situtated at the center of a particle. Otherwise every cell with a "phaseFractionField" smaller than one will be considered.

* "direction": requires an integer value. If set to 0, samples will be drawn in the x-direction. If set to 1, samples will be drawn in the y-direction. If set to 2, samples will be drawn in the z-direction.

NOTE: if "save2Bin" is true, two contiguous binning operations are required! The first one will bin the numerator and the second the denominator.
NOTE: if "save2Disk" is true two files will be dumped. The first containing numerators and the second denominators.


Example
-------
```
...
                 
   "2PC": {
                  "type": "twoPointsCorr",
                  "fieldToSample"       : "U",
                  "filteredField"       : "U_Favre",
                  "direction"           : 0,
                  "phaseFractionField"  : "interFace_Favre",
                  "sampleCount"         : 120,
                  "sampleDelta"         : 0.1,
                  "save2Disk"           : true,
                  "save2Bin"            : true,
                  "lagrangian"          : false
              },
...
```

Go back
-----------
 - [sampling](12_sampling.md)
