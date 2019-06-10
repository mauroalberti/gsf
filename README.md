
# gsf 
**gsf** is a library for the processing of structural geology data. It is developed in Python 3 (module *pygsf*).

## Summary

Currently the module allows to create several geometric and geological "objects" and to process them.

Geometric classes comprise cartesian points (**Point**), vectors (**Vect**) and planes (**CPlane**).

Geological direction and axes (**Direct** and **Axis**), as well geological planes (**CPlane**), are expressed via usual geological notation, i.e., dip direction and dip angle.

For the processing of fault data, it is possible to create slickenline-type objects (**Slick**), with or without known movement sense. 

A slickenline can be combined to a geological fault to create a fault - slickenline datum (**Fault**), from which it is possible to derive P, T and B axes, and the M plane (**PTBAxes**). 

A quaternion module (**Quaternion** class) allows to rotate geological data.

## Installation as a Python 3 module

You should have *python-setuptools* already installed, if not you can install in a Linux Mint shell via apt or pip3, e.g.:
```
sudo apt-get install python-setuptools
```
To install *gsf*, from the shell you can run in Linux Mint:
```
sudo python3 setup.py install
```
You can then check if the module is installed in Python3 by successfully importing the module via:
```python
import pygsf
```
## Docs

There are a few Jupyter notebooks describing some of the features of the module in the docs/notebooks folder.b)





