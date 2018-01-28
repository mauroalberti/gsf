# -*- coding: utf-8 -*-


from ..ptbaxes import *


"""

Dataset from: 
Kagan Y.Y., 1991
3-D rotation of double-couple earthquake sources
 
   values in tuples are (trend, plunge)
   
"""

# original values

k91_fs_T_axis_vals, k91_fs_P_axis_vals = (120, 24), (232, 41)
k91_ss_T_axis_vals, k91_ss_P_axis_vals = (295, 55), ( 51, 17)


# gsf-scope derived parameters

k91_fs_T_gaxis, k91_fs_P_gaxis = GAxis(*k91_fs_T_axis_vals), GAxis(*k91_fs_P_axis_vals)
k91_ss_T_gaxis, k91_ss_P_gaxis = GAxis(*k91_ss_T_axis_vals), GAxis(*k91_ss_P_axis_vals)

k91_fs_PTBaxes = PTBAxes(
    p_axis=k91_fs_P_gaxis,
    t_axis=k91_fs_T_gaxis)

k91_ss_PTBaxes = PTBAxes(
    p_axis=k91_ss_P_gaxis,
    t_axis=k91_ss_T_gaxis)

k91_fs_quater = k91_fs_PTBaxes.to_quaternion()
k91_ss_quater = k91_ss_PTBaxes.to_quaternion()


