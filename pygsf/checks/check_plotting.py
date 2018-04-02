
from pygsf.geometry import *
from pygsf.geotype_checks import *
from pygsf.plotting.stereonets import splot

gv1 = GVect(312, -45)
gv2 = GVect(300, -20)
gv3 = GVect(180, 10)

ga1 = GAxis(180, 10)


splot([(gv1, gv2, "m=s, c=red"), gv3, ], force='lower')

