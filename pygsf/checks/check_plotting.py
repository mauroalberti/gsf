
from pygsf.geometry import *
from pygsf.geotype_checks import *
from pygsf.plotting.stereonets import splot

ga1 = GAxis(312, -45)
ga2 = GAxis(300, -20)
ga3 = GAxis(180, 10)

splot([(ga1, "m=s,c=blue"), ga2, ga3], force='')

