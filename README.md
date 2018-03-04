
[![Build Status](https://travis-ci.org/mauroalberti/gsf.svg?branch=master)](https://travis-ci.org/mauroalberti/gsf)

# gsf 
**gsf** is a library for the processing of structural geology data. It is developed in Python 3 (module *gsf_py*), and also as a preliminary and in-progress Haskell version (module *gsf_hs*).

## Summary

Currently the module allows to create several geometric and geological "objects" and to process them.

Geometric classes comprise cartesian points (**Point**), vectors (**Vect**) and planes (**Plane**).

Geological direction and axes (**GVect** and **GAxis**), as well geological planes (**GPlane**), are expressed via usual geological notation, i.e., dip direction and dip angle.

For the processing of fault data, it is possible to create slickenline-type objects (**Slickenline**), with or without known movement sense. 

A slickenline can be combined to a geological fault to create a fault - slickenline datum (**FaultSlick**), from which it is possible to derive P, T and B axes, and the M plane (**PTBAxes**). 

A quaternion module (**Quaternion** class) allows to rotate geological data.

## Installation as a Python 3 module

You should have *python-setuptools* already installed, if not you can install in a Linux Mint shell via apt or pip3, e.g.:
```
sudo apt-get install -y python-setuptools
```
To install *gsf*, from the shell you can run in Linux Mint:
```
sudo python3 setup.py install
```
You can then check if the module is installed in Python3 by successfully importing the module via:
```python
import gsf_py
```

## Usage examples for the Python version 

*(In progress, mirrors the content of the Jupyter notebook in the folder "Notebook")*

Since geological data can be expressed in terms of geometric and geographical parameters,
the basics of this module are geometric concepts: points, planes and vectors. From these concepts, specialised geological concepts are derived: geological vectors and axes, and geological planes. From the fundamental side, geological vectors and axes are vectors, while geological planes are geometric planes. The difference in more in the way to express the values of these structures. In structural geology, orientations are expressed via angles from references directions (polar coordinates), such as the North or as the local horizontal plane. In the geometric realm, orientations are mainly expressed as Cartesian coordinates.
We start by considering Cartesian points.

### Cartesian points and planes 

A point can be created in the usual way:


```python
from gsf_py.geometry import *
```


```python
p1 = Point(1.0, 2.4, 0.2)  # definition of a Point instance
```

Among other properties, we can calculate its distance from the reference frame origin via the method abs(): 


```python
abs(p1)  # distance of a point from the origin
```




    2.6076809620810595



Given another point, we can calculate the 3D and the horizontal distance (2D) between two points.


```python
p2 = Point(0.9, 4.2, 10.5)
```


```python
p1.dist_3d(p2)  # 3D distance between two points
```




    10.45657687773585




```python
p1.dist_2d(p2)  # horizontal (2D) distance between two points
```




    1.8027756377319948



Other possibilities are to translate the point via a triad of cartesian values or directly via a vector, to check if two points are within a given range and to convert a point to a vector.

A Cartesian plane can be defined in a few ways:
    from three points:


```python
pl1 = Plane.from_points(Point(0, 0, 0), Point(1, 0, 0), Point(0, 1, 0))  # definition of a plane from three points
```


```python
print(pl1)
```

    Plane(0.0000, 0.0000, 1.0000, 0.0000)


Those returned are the four coefficient (a, b, c and d) defining the Cartesian plane, defined by the equation: 

   *ax + by + cz = d*

It can be seen that for the provided example the equation is satisfied for all *x* and *y* values when *z* is always zero, i.e. the Cartesian plane is a horizontal plane passing through the frame origin.
    

We calculate the versor normal to this plane:


```python
normal_versor = pl1.nversor()  # versor (unit vector) normal to the provided Cartesian plane
```


```python
print(normal_versor)
```

    Vect(0.0000, 0.0000, 1.0000)


And we see that the versor is the vertical axis, as expected.

As another example, we calculate the intersection between two vectors, expressed by a versor.


```python
pl1, pl2 = Plane(1, 0, 0, 0), Plane(0, 0, 1, 0)
inters_v = pl1.inters_versor(pl2)  # intersection versor between two Cartesian planes 
print(inters_v)
```

    Vect(0.0000, -1.0000, 0.0000)


### Vector

Vector creation and manipulation are straightforward:


```python
from gsf_py.geometry import *
v1, v2 = Vect(3.1, 7.2, 5.6), Vect(4.2, 9.17, 8.0)
```


```python
v1 + v2  # vector addition
```




    Vect(7.3000, 16.3700, 13.6000)




```python
v1 - v2  # vector subtraction
```




    Vect(-1.1000, -1.9700, -2.4000)



Scalar and vector products are obtained via:


```python
v1.sp(v2)  # scalar product
```




    123.84399999999999




```python
v1.vp(v2)  # vector product
```




    Vect(6.2480, -1.2800, -1.8130)



The angle (in degrees) between two vectors, and the check if they are sub-parallel or sub-orthogonal: 


```python
v1.angle(v2)  # angle in degrees bwtween two Cartesian vectors
```




    3.0646173501805807




```python
v1.almost_parallel(v2)  # is v1 sub-parallel to v2?
```




    False




```python
v1.is_suborthogonal(v2)  # is v1 sub-orthogonal to v2?
```




    False



A vector can be converted to a geological vector or axis by using the gvect() and gaxis() methods.


```python
gv1 = v1.gvect()  # conversion from Cartesian vector to geological vector
```


```python
print(gv1)
```

    GVect(023.29, -35.54)



```python
ga2 = v2.gaxis()  # conversion from Cartesian vector to geological axis
```


```python
print(ga2)
```

    GAxis(024.61, -38.42)


### Geological vectors, axes and planes

A *geological vector* is a vector in the 3D space with unit length and a direction defined by a trend (from the North, 0°-360°) and a plunge (-90° to 90°, where positive values are downward-directed while negative ones are upward-directed).


```python
gv1, gv2 = GVect(312, 45), GVect(92, -38)  # gv1 and gv2 are two geological vectors defined by trend and plunge values
```

*Geological axes* are similar to geological vectors, but do not have a specific direction, i.e., they have only an orientation. As for geological vectors, they are defined by a trend and a plunge, but the two possible, opposite directions are both allowed and considered for calculations (e.g., the angles between axes), even if they are not explicit.

We can create a geological axis given trend and plunge values, or convert from geological vector to a geological axis: 


```python
ga1 = GAxis(219, 24)  # creating a geological axis given trend and plunge
```


```python
ga2 = gv2.as_axis()  # converting a geological vector to a geological axis
```


```python
print(ga1, ga2)
```

    GAxis(219.00, +24.00) GAxis(092.00, -38.00)


A *geologic plane* is a plane with orientations expressed via geological convention, i.e. the azimuth from North of the strike or dip direction, and the dip angle. The module follows the dip direction convention:


```python
gpl1 = GPlane(112, 67)  # definition of a geological plane instance
```


```python
print(gpl1)
```

    GPlane(112.00, +67.00)


As previously said, the main distinction between geological axes and vectors is that the the former have not a direction, while the latter are oriented. This difference is reflected for instance in the calculation of the angle between two vectors and two axes:


```python
vector_angle = gv1.angle(gv2)  # angle (in degrees) between two geological vectors 
```


```python
print(vector_angle)
```

    149.56272807646775


We convert the geological vectors to axes:


```python
ga1, ga2 = gv1.as_axis(), gv2.as_axis()
```


```python
axis_angle = ga1.angle(ga2)  # angle (in degrees) between two geological axes
```


```python
print(axis_angle)
```

    30.43727192353225


The angle between the two axes is the complement to 180° of the angle between the two geological vectors.

In addition to opposite, downward and upward geological vectors, it is possible to calculate the geological plane common to two geological vectors:


```python
gplane = gv1.common_plane(gv2)  # geological plane common to two geological vectors (gv1 and gv2)
```


```python
print(gplane)
```

    GPlane(310.64, +45.01)


where the former value is the dip direction of the plane and the latter is the dip angle (downward since positive)

as well as the vector normal to both geological vectors:


```python
ngv = gv1.normal_gvect(gv2)  # geological vector normal to gv1 and gv2 geological vectors
```


```python
print(ngv)
```

    GVect(130.64, +44.99)


where the first value is the dip direction and the second one is the dip angle.

Considering just a single geological vector, the geological plane normal to the vector is obtained via:


```python
ngp = gv1.normal_gplane()  # geological plane normal to a given geological vector
```


```python
print(ngp)
```

    GPlane(132.00, +45.00)

