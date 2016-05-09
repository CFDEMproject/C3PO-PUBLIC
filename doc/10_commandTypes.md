Command types
==============

_Command types_ are used to identify the specific container for each command.
There can be three _command types_:

* _operation_

* _selector_

* _control_

Command type `operation`
===============

Syntax
---------------------
Commands of type `operation` have to be followed by one of the following _command specifiers_:

* _filtering_ 

* _sampling_

* _binning_

They are separated using spaces. Multiple _operation command types_ can be used in a _c3po.input_ file.

Description
---------------------
The _operation command type_ tells CPPPO to perform an operation on the data. The subsequent _command specifier_ declares which operation has to be performed.

Example
---------------------
In `c3po.input`:

```
...

operation filtering filtU
operation filtering filtP

operation sampling uSlip
operation sampling betaP
operation binning  uSlip
operation binning  betaP

...

```
These lines will run six different operations.

Command type `selector`
===============

Syntax
---------------------
The _selector command type_ needs to be specified using one of the following commands:

* _cellIJK_ 

* _cellUnstruct_ 

* _bubble_

* _filter_

They are separated using spaces. Multiple _selector command types_ can be used in a _c3po.input_ file.

Description
---------------------
The _selector command type_ tells CPPPO how to calculate the number of cells or particles within the filter, i.e., how to select the required cells or particles. In order to separate the specification of the size of the region and the algorithm to select the cells in this region, a pair of _selector command type_ declarations needs to be used: (i) the `filter` specifier is used to define the filter size, while (ii) the subsequent _command specifier_ (either `cellIJK` or `particleANN`) declares which algorithm to select the cells or particles has to be used.

Example
---------------------
In `c3po.input`:

```
...

selector  cellIJK       myCell
selector  bubble        myBubble

selector  filter        myFilter0
selector  filter        myFilter1

...

```
These lines will activate two pairs of selectors. As explained above, this is because a _selector command type_ with a _filter command specifier_ needs a corresponding _selector command type_ with a _cellIJK_ or _particleANN command specifier_ to run. 

Control  command type
===============

Syntax
---------------------
The _control command type_ needs to be specified using one of the following commands:

* _run_ 

They are separated using spaces. Only one _control_ command is allowed in a cpppo.input file.

Description
---------------------
The _control command type_ command will make CPPPO run as a stand-alone (i.e. without any interface) for a certain time.
The run time must be specified after _run_ instead of the _command name_. 

Example
---------------------
In `c3po.input`:

```
...

control   run 2
...

```
This line will run CPPPO as standalone for 2 seconds.

Go back
-----------
 - [main](01_main.md) 

