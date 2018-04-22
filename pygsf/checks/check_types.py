

from pygsf.geography import *

gv1 = GVect(100, 23)

print(isinstance(gv1, GVect))
print(isinstance(gv1, GAxis))


ga1 = GAxis(100, 23)

print(isinstance(ga1, GVect))
print(isinstance(ga1, GAxis))
