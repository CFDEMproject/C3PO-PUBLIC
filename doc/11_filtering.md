Filtering command specifiers
======================


_Filtering command specifiers_always need to be followed by the _command name_. Any _operation_ declared using a _filtering command specifier_ is called _filtering operation_ and has to be defined in the c3po.json file.
The _filtering command specifier_ refers to the functional form of the filter Kernel function. At the moment, only the Top-Hat is available:

* _filteringTopHat_

Syntax  
-------
```
...

"favreU": 
            {
              "mainSettings":
              {
               "VectorfieldsToFilter": "U ",
               "ScalarfieldsToFilter": " ",
               "lagrangian": false             
              },
              
              "varianceSettings":
              {
               "VectorfieldsForVarianceName1": "U",
               "VectorfieldsForVarianceComputeOffDiagonal": [false],
               "ScalarfieldsForVectorScalarMixedVariance"  : " off",
               "ScalarfieldsForVarianceName1" : "",
               "ScalarfieldsForVarianceName2" : "",
               "computeVariance":true
             
              },       
                    
              "kernelSettings":
               {
                "weightFields":"",
                "inverseWeightFields":"alpha"
               }
            },
...

```
Every _filtering operation_ requires the user to fill some fields in the c3po.json file according to the "type" entry. However, any _filtering operation_ requires the following fields:

* "mainSettings":     this sub-dictionary requires several entries:
 * "VectorfieldsToFilter": requires a string of characters. It has to be filled with the names of the vector fields to filter. 
 * "ScalarfieldsToFilter": requires a string of characters. It has to be filled with the names of the scalar fields to filter. 
 * "lagrangian":           requires a boolean value. If set to _false_, every cell in the domain will be filtered, otherwise only cells at particle/probe centres will be filtered (Note that the "lagrangian": false algorithm is faster). If this is set to _true_ the name of the probes/particles group should be provided with an additional entry:
 * "probesName": requires a string of characters. The user has to enter the probes/particles group that will be used by this operation. Note that just one group per operation can be defined.

* "kernelSettings":  this sub-dictionary allows to customize the Kernel function adding an arbitrary number of weights. It requires several entries:
 * "weightFields" : requires a string of characters. It has to be filled with the names of the scalar fields to use as weights.
 * "inverseWeightFields": requires a string of characters. It has to be filled with the names of the scalar fields to use as weights. The weight value will be one minus the field value


Additionally, filtering operation allows to compute variance fields. In order to do that, several new entries need to be provided:

* "varianceSettings": this sub-dictionary requires several entries:
 * "VectorfieldsForVarianceName1" : requires a string of characters. This is the list of vector fields for variance calculation.
 * "ScalarfieldsForVarianceName1" : requires a string of characters. This is the list of scalar fields for variance calculation.
 * "VectorfieldsForVarianceComputeOffDiagonal" : requires an array of bools. If any element set to 'false', only the diagonal components of the variance field for the corresponding field in  "VectorfieldsForVarianceName1" will be computed.
 * "ScalarfieldsForVectorScalarMixedVariance" : requires a string of characters. This is the list of scalar fields for vector-scalar correlation. Enter "off" to disable.
 * "VectorfieldsForVarianceName2" : requires a string of characters. This is the list of vector fields for vector-vector correlation.
 * "ScalarfieldsForVarianceName2" : requires a string of characters. This is the list of scalar fields for scalar-scalar correlation.
 * "computeVariance": requires a bool. If set to false CPPPO will skip the variance calculation.

Example  
-------
In `c3po.input`

```
...

operation filteringTopHat favreU

...

```

In `c3po.json`

```
...

"favreU": 
            {
              "mainSettings":
              {
               "VectorfieldsToFilter": "U ",
               "ScalarfieldsToFilter": " ",
               "lagrangian": false             
              },
              
              "varianceSettings":
              {
               "VectorfieldsForVarianceName1": "U",
               "VectorfieldsForVarianceComputeOffDiagonal": [false],
               "ScalarfieldsForVectorScalarMixedVariance"  : " off",
               "ScalarfieldsForVarianceName1" : "",
               "ScalarfieldsForVarianceName2" : "",
               "computeVariance":true
             
              },       
                    
              "kernelSettings":
               {
                "weightFields":"",
                "inverseWeightFields":"alpha"
               }
            },
...

```


Go back
-----------
 - [main](01_main.md) 

