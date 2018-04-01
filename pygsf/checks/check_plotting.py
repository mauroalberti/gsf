
from pygsf.geometry import *
from pygsf.geotype_checks import *
from pygsf.plotting.stereonets import splot

gv1 = GVect(312, -45)
gv2 = GVect(90, 10)

ga1 = GAxis(180, 10)

data = [gv1, gv2, ga1]

gvects = list(filter(is_gvect, data))
gaxes = list(filter(is_gaxis, data))

gvects_notup = list(filter(is_not_upward, gvects))
gvects_up = list(filter(is_upward, gvects))
gaxes_notup = list(filter(is_not_upward, gaxes))
gaxes_up = list(filter(is_upward, gaxes))



splot([gv1, gv2, ga1])

