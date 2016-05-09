CPPPO input files
===============


Description
---------------------
CPPPO needs no coding from the user. Instead, CPPPO processes text files containing _commands_ and basic mesh information, which have to be provided by the user. Every _command_ handed over to CPPPO (as specified in the main input scipt) needs to be accompanied with a corresponding file contained in the folder `./c3po_control`. The three key input files to CPPPO are:

- `c3po.input`
- `c3po.json`
- `mesh.json`


Each _command_ is associated to a _command type_, a _command specifier_ and a _command name_ with the exception of the _control run_ command, which does not require any _command name_. _command_ declarations are made in the _CPPPO.input_ file, while the accompanying definitions are specified in the `c3po.json` file. Differences between declaration and definition are that:

* Declaration: they are used to inform CPPPO about which and how many objects have to be created, as well as their _names_, i.e., a word used to identify a particular _command_. The name of each _command_ will be used also to identify CPPPO's output (e.g., the file name). The declarations can be used to activate or deactive certain commands (i.e., just comment-out a specific declaration using a `#`).

* Definition: once every _command_ is declared and objects are created in CPPPO core, user defined parameters need to be provided. This is done via the definitions provided in the `c3po.json` file.

`c3po.input`
---------------------
CPPPO reads _command types_, _command specifiers_ and _command names_ from the _c3po.input_ file. There are three different _command types_. 

1.  _operation_

2.  _selector_

3.  _control_

###Example

Each type can be used multiple times with different _command specifiers_ and/or _command names_.

```
operation    filteringTopHat        favreU
operation    filteringTopHat        algebraicP

operation    samplingGeneral        2pointCorr

operation    binning                2PointCorrBins

selector     filter                 myFilter0
selector     cellIJK                cellSelector0

selector     filter                 myFilter1
selector     cellIJK                cellSelector1   


```
_command types_, _command specifiers_ and _command names_ are separated using at least one space character. 
NOTE: _commands_ are read and numbered sequentially in CPPPO core. This means that the user has to be sure that if two different _commands_ need to interact, they are addressed with the same id by CPPPO core. In the above example, _operation sampling 2pointCorr_ and _operation binning 2PointCorr_ interact (they are named with the same _command name_ in order to make it clear) and thus, they are both the first _operations_, among those with the same _command specifier_, declared in _CPPPO.input_. This is the only ordering required when writing a _CPPPO.input_ file. Thus, the example below produces the exact same behavior as the previous one:


```
operation    filteringTopHat        favreU
operation    filteringTopHat        algebraicP

operation    samplingGeneral        2pointCorr

operation    binning                2PointCorrBins

selector     cellIJK                cellSelector0
selector     cellIJK                cellSelector1   

selector     filter                 myFilter0
selector     filter                 myFilter1

```
`c3po.json`
---------------------
The _c3po.json_ file is used to define the _commands_ handed over to CPPPO, and to define other basic settings for CPPPO. The "mainSetting" field provides the latter functionality.

###Example

```
{
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
        "Scalarfields": " p alpha "
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
  
   "myFilter0":	  { 
                     "CoordSys": 0, 
                     "x":  1.0 , 
                     "y":  1.0 ,
                     "z":  1.0 
                  },
    "myFilter1": { 
                     "CoordSys": 0, 
                     "x":  2.0, 
                     "y":  2.0,
                     "z":  2.0 
                  },
    
   "favreU": 
            {
             "mainSettings":
              {
               "VectorfieldsToFilter": "U ",
               "ScalarfieldsToFilter": " ",
               "lagrangian": false             
              },
                    
              "kernelSettings":
               {
                "weightFields":"",
                "inverseWeightFields":"alpha"
               }
            },
   
   "algebraicP": 
             {
              "mainSettings":
              {
               "VectorfieldsToFilter": "U ",
               "ScalarfieldsToFilter": " ",
               "lagrangian": false
             
              },
                    
              "kernelSettings":
               {
                "weightFields":"",
                "inverseWeightFields":""
               }
            }

            
...

}

```
The above example shows a part of the _c3po.json_ file related to the previous _c3po.input_ file. As can be seen, every _command_ is defined in the _c3po.json_ file using its _command name_ declared in the _c3po.input_ script. In addition, the "mainSettings" field has to be specified, and has to include the following subfields:

* "doFiltering": requires a boolean value. If set to false, CPPPO will skip any _filtering operation_. 
* "doSamling":   requires a boolean value. If set to false, CPPPO will skip any _sampling operation_. 
* "doBinning":   requires a boolean value. If set to false, CPPPO will skip any _binning operation_.
* "dumpFormat":  requires a string of characters. Allows to set the default output format for CPPPO and it can be set either to "json" or "hdf5". Note, "hdf5" only works if CPPPO is compiled with the HDF5 option.
* "verbose": requires a boolean value. If set to true, CPPPO will dump information regarding performed operations to screen.
* "FieldsToRegister": requires a string of characters for every sub-entry. Here the user has to specify which fields, among those available from the interface (e.g., the simulator or the .csv file), have to be registered in CPPPO core. There is also a distiction between vector and scalar fields.
* "filterBCs" : requires a string of characters for every sub-entry. Here the user has to specify the boundary conditions to use for filtering operations. Filtering boundary conditions are symmetric and have to be specified for every direction in a cartesian reference frame. The user can enter "periodic" or "wall".
* "storageWriteFields" and "storageWriteParticles" : switch to activate dumping of C3PO's main storage to disk.
* "storageWriteFields" and "storageWriteParticles" : switch to activate dumping of C3PO's interface to disk. If you are using the CSV interface, CPPPO will dump filtered fields in CSV format while, if you are using the OpenFOAM interface, CPPPO will dump filtered fields in OpenFOAM format (as any other field).

Some specific functions may require additional entries. In this case, the user should refer to the corresponding documentation page.

`mesh.json`
---------------------
###Example

```
{
    "mainSettings": 
     {
        "checkTolerance":   1e-07,
        "filterWidthTolerance": 1e-03,
        "verbose":          false
        
     }



}
```
The `mesh.json` file contains information regarding the computational grid, as well as the checks that are performed on this grid. The number of information needed varies from OpenFOAM applications to .csv application. The user shoud supply the following information in order to control mesh checks:

* "checkTolerance":         this controls the relative tolerance for the cell volume checks (recommended value: 1e-5)
* "filterWidthTolerance":   this controls the absolute tolerance (in SI units) for detecting whether a cell is within the filter or not. The optimum value varies from case to case and should be a fraction of the minimum cell legth.

However, one field is always required:

* "verbose": requires a boolean value. If set to true CPPPO will dump information regarding mesh and _selecting operations_ to screen.

In addition, more entries are required when using the CSV interface.

 
Go back
-----------
 - [main](01_main.md)

