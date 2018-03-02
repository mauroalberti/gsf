
[![Build Status](https://travis-ci.org/mauroalberti/gsf.svg?branch=master)](https://travis-ci.org/mauroalberti/gsf)

# gsf 
**gsf** is a library for the processing of structural geology data.

It is developed in Python 3, and a Haskell version is also in development.

### Summary

Currently the module allows to create several geometric and geological "objects" and to process them.

Geometric classes comprise cartesian points (**Point**), vectors (**Vect**) and planes (**Plane**).

Geological direction and axes (**GVect** and **GAxis**), as well geological planes (**GPlane**), are expressed via usual geological notation, i.e., dip direction and dip angle.

For the processing of fault data, it is possible to create slickenline-type objects (**Slickenline**), with or without known movement sense. 

A slickenline can be combined to a geological fault to create a fault - slickenline datum (**FaultSlick**), from which it is possible to derive P, T and B axes, and the M plane (**PTBAxes**). 

A quaternion module (**Quaternion** class) allows to rotate geological data.

