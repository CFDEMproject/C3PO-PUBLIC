CPPPO post-filtering operations
===============


Description
---------------------
Post-filtering operations are those operations performed on filtered fields which generate new fields (i.e. computation of gradients or shear rates).
At the current state, gradient and shear rate computation is available.

Syntax
-----------

In order to calculate gradients of the filtered fields, it is necessary to specify an additional entry in 'c3po.json' inside the "fieldsToRegister" object in "mainSettings":
* "computeScalarGradient" : requires a string of characters. These are the name of the filtered scalar fields to process. The output will be in the form gradFieldName Note _the name of the filtered field is required_, which is in the form: originalFieldName\_operationName\_filterName.
* "computeVectorGradient" : requires a string of characters. These are the name of the filtered vector fields to process. The output will be in the form gradFieldNameComponent and it consists of three fields each one referreng to the _gradient of a particular component_ of the vector field. Note _the name of the filtered field is required_, which is in the form: originalFieldName\_operationName\_filterName.

* "computeShearRate" : requires a string of characters. These are the name of the filtered vector fields to process. The output will be in the form shearFieldNameComponent. Note _the name of the filtered field is required_, which is in the form: originalFieldName\_operationName\_filterName.

Example
-----------
In this example also second derivatives are computed.
```
"mainSettings":
   {
        "doFiltering": true,
        "doSampling":  true,
        "doBinning":   true,
        "dumpFormat": "hdf5",
        "verbose": true,
        "FieldsToRegister":
        {
        "Vectorfields": " U ",
        "Scalarfields": " p alpha ",
        "computeScalarGradient" : "p_favre1_myFilter0 p_favre1_myFilter1",
        "computeScalarGradient" : "U_favre1_myFilter0 U_favre1_myFilter1 gradU_favre1_myFilter0X",
        "computeSchearRate" : "U_favre1_myFilter0 U_favre1_myFilter1"       
        },
        "filterBCs":
        {
        "x":"periodic",
        "y":"periodic",
        "z":"periodic"
        },
        "storageWriteFields"    :  false,
        "storageWriteParticles" :  false,
        "interfaceWriteFields" :  false
     
   },
```

Go back
-----------
 - [input](02_c3poInput.md)

