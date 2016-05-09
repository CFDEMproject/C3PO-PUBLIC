Sampling formula
======================


Description
----------------------

CPPPO allows the user to perform basic algebraic operations on sampled fields and to register the result with the operation marker. A _formula_ is an optional entry for every _sampling operation_.

* "Formula": requires a string of characters. The user can type the formula she/he wants CPPPO to execute.

The following rules applies to user defined formulas:

* The user can specify a numerator and a denominator with multiple operations for each. The expression for the numerator and denominator MUST be in brackets and separated with the '%' symbol.
* CPPPO will fill the data computed via the formular into a sample field with the extension 'formula0'
* elements specified via the keyword 'VFieldsToSample' can be accessed with 'vec1', 'vec2', ..., 'vec9'. It is not possible to use more than 9 vectors (or scalars) in such a way
* elements specified via the keyword 'SFieldsToSample' can be accessed with 'scalar1', 'scalar2', ..., 'scalar9'.
* Operations in the numerator and the denominator are executed from left to right. Thus, the order of the operations defines the output, and brackets in the numerator and denominator will be NOT respected!
* Thus, the formula "(vec1 - vec2) * scalar1" will give a different result as "scalar1 * (vec1 - vec2)"  
* The keyword "saveOnlyFormula" can be specified in case a formula is specified. If set to true, only the result of the formular will be saved (to bin or disk), but NOT the sampled vector and scalar fields.
* The user can also directly refer to sampled fields using numbers from 0 to 9. CPPPO will order the fields considering vector fields first using the ordering provided by the user input.
* The only operations allowed are:
 * '*' : multiplication
 * '/' : division
 * '+' : sum
 * '-' : difference

Note that the computed value will be considered as a new sample and, if "save2Bin" is set to 'true', it will require an additional _binning operation_. The user may want to set "saveOnlyFormula" : true in order to avoid specifying an uncessary large amount of bins.

The formula class also allows the user to normalize each sampled value using any normalization function provided in 'core/multiphase_FlowBasic.H'. At the current state, only
interphase exchange coefficients are availble as normalizing functions. Additional entries required are:
* "NormalizeFormula" :  this is the Json object that must contain the following entries:
 *  "Normalization" : requires a string of characters. The user must specify the kind of normalization as well as the required fields. The first character should be the normailzation Id, while the remaining ones are the required fields. These are specified using 'm', 's' or 'v' is they are markers, scalar of vector fields which _have been already registered in the sampling operation_ followed by a number indicating the order of declaration. Available normalization Ids are:
  * D : is calling for the interphase exchange coefficient of a drag law. It requires further fields:
   * "dp" : requires a double. It is the particle diameter.
   * "rhoP" : requires a double. It is the particle density.
   * "g" : requires a double. It is the gravitational acceleration.
   * "etaFluid": requires a double. It is the fluid kinematic viscosity.
   * "rhoFluid" : requires a double. It is the fluid density.
   * "dragLaw" : requires an integer. It is the drag law Id ( 0 = Beetstra , 1 = WenYu , 2 = KochHill  ). Every drag law requires three fields to be provided in the following order: phase1 fraction, phase2 velocity, phase1 velocity. 

The formula class allows to introduce limits to the sampled values in order to get better statistics. The use of limiters requires new entries:

* "formulaLimiterNumerator" : requires an array of doubles. The first element is the minimum and the last one is the maximum for the numerator.
* "formulaLimiterDenominator" : requires an array of doubles. The first element is the minimum and the last one is the maximum for the denominator.
* "formulaLimiterSymmetry" :  requires an array of bools. The first element relates to the numerator and the second to the denominator. If set to true, the limiter would be considered in terms of absolute values.

The formula class allows the absolute values of the numerator, denominator or both to be used in the formula. To sample absolute values, the following entries are required:

* "absoluteNumerator" : requires an bool. If set to true, the absolute value of the numerator will be used in the formula. Default value is false
* "absoluteDenominator" : requires an bool. If set to true, the absolute value of the denominator will be used in the formula. Default value is false 

Example
-------
Sampling operation with formula to sample the drag coefficient.
```
...
                 
   "sample0": {
                  "type"            : "general",
                  "marker"          : "solidsFraction_filtFavre",
                  "VFieldsToSample" : "fluidVel_x_filtFavreGas solidVel_x_filtFavre force_x_filtFavre",
                  "Formula"        : "( vec3 ) % ( vec1 - vec2 ) ",
		  "absoluteNumerator": true,	
		  "absoluteDenominator": true, 
                  "formulaLimiterNumerator"   : [0,  1e99],
                  "formulaLimiterDenominator" : [0.1,  1e99],
                  "formulaLimiterSymmetry"    : [true, true],
                  "NormalizeFormula"          : 
                                              {

                   "Normalization" : "Dm1v1v2",
                   "dp"   : 75e-6,
                   "rhoP" : 1000, 
                   "g"    : 9.81, 
                   "etaFluid" :     1.8e-5, 
                   "rhoFluid" :     1.3, 
                   "dragLaw"  :     1 
                  },

                  "saveOnlyFormula" : true,
                  "component"       : 0,
                  "sampleCount"     : -1,
                  "save2Disk"       : true,
                  "save2Bin"        : true,
                  "lagrangian"      : false,
                  "overwrite"       : true                  
               },

...
```
In this example the formula entered corresponds to:
```
 (force_x_filtFavre) / (fluidVel_x_filtFavreGas - solidVel_x_filtFavre)  
```

Go back
-----------
 - [sampling](12_sampling.md)
