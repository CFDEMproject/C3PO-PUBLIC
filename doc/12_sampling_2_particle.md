Operationsampling particles
======================


"type": "particles"
----------------------
"Particles" _sampling operations_ are used to extract additional quantities from a lagrangian field or filtered quantities at locations
of interest. The _command specifier_ is:

* _samplingParticles_

This sampling operation requires "lagrangian":true.
If "sampleCount" is equal to -1, all the particles are sampled (they are still limited by "selective").

Syntax
-----
```
...

 "parsampleS": {
                  "marker"          : "Reynolds SauterW",
                  "VFieldsToSample" : "U_favre ",
                  "component"       : 0,
                  "SFieldsToSample" : "particleForce Sherwood Reynolds ",
                  "maxParticleRadius":0.6,
                  "minParticleRadius":0.5,            
                  "filteredVelocityForRe":"U_favre",
                  "nu":0.1,
                  "Pr":1,
                  "filteredXForSh":"Y_favre",
                  "Xref":2.0,
                  "Qid":0,
                  "particleFractionTot":"inside_average",
                  "particleFraction_i":["fractionS_average","fractionS_average"],
                  "particleDiameter_i":[1.0,2.0],
                  "fraction_id": 0,
                  "sampleCount"     : -1,
                  "save2Disk"       : true,
                  "save2Bin"        : true,
                  "lagrangian"      : true,
                  "probesName" :"polyPar",
                   "selective":true,
                   "max":[10.5,  8,  8],
                   "min":[ 2.5,    0,  0], 
                  "overwrite"       : false                        
               },


...

This sampling operation allows to compute quantities of interest in the field of multiphase flows.
The following entries can be used as both markers or samples:

* "Reynolds" : this computes the particle Reynolds number. It requires the user to provide additional entries:
 * "filteredVelocityForRe" : requires a string of characters. It is the name of the velocity field used to compute the Reynolds number. The field magnitude is used.
 * "nu" : requires a double. It is the fluid kinematic viscosity.
 * "particleFractionTot" : requires a string of characters. It is the name of the total particle fraction to be used.

* "Sherwood" : this computes the particle Sherwood number. It requires the user to provide additional entries:
 * All the entries required for "Reynolds".
 * "Pr" : requries a double. It is the value of the Prandtl number.
 * "filteredXForSh" : requires a string of characters. It is the name of the filtered scalar field used to compute the Sherwood number.
 * "Xref" : requires a double. It is the particle surface scalar value.
 * "Qid" : requires an integer. This number refers to the particle scalar representing the total dimensionless fluid-particle heat transferred. 
 
* "SauterW" : this computes the particle Sauter diameter. It requires the user to provide additional entries:
 * "particleDiameter_i" : requires an array of doubles. It is the list of particle diameters.
 * "particleFraction_i" : requires an array of strings. It is the list of names of the particle fraction fields corresponding to the "particleDiameter_i".
 * "fraction_id" : requires an integer. Is the id of the current particle class according to "particleDiameter_i".
 * "particleFractionTot" : requires a string of characters. It is the name of the total particle fraction to be used.
 
* "particleForce" : returns the magnitude of the total force acting on the particle.

In addition the user can specify:

* "maxParticleRadius" : requires a double. Is the maximum particle radius to consider.
* "minParticleRadius" : requires a double. Is the minimum particle radius to consider.

```
Example
-------
In `c3po.input`:
```
...

operation samplingParticles  parsampleS

...
```
In `c3po.json`:
```
...
                 
  "parsampleS": {
                  "marker"          : "Reynolds SauterW",
                  "VFieldsToSample" : "U_favre ",
                  "component"       : 0,
                  "SFieldsToSample" : "particleForce Sherwood Reynolds ",
                  "maxParticleRadius":0.6,
                  "minParticleRadius":0.5,            
                  "filteredVelocityForRe":"U_favre",
                  "nu":0.1,
                  "Pr":1,
                  "filteredXForSh":"Y_favre",
                  "Xref":2.0,
                  "Qid":0,
                  "particleFractionTot":"inside_average",
                  "particleFraction_i":["fractionS_average","fractionS_average"],
                  "particleDiameter_i":[1.0,2.0],
                  "fraction_id": 0,
                  "sampleCount"     : -1,
                  "save2Disk"       : true,
                  "save2Bin"        : true,
                  "lagrangian"      : true,
                  "probesName" :"polyPar",
                   "selective":true,
                   "max":[10.5,  8,  8],
                   "min":[ 2.5,    0,  0], 
                  "overwrite"       : false                        
               },
...
```

Go back
-----------
 - [sampling](12_sampling.md)
