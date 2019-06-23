# -*- coding: utf-8 -*-


from typing import List, Dict, Union

import numpy as np


def get_statistics(vals: Union[List, np.array]) -> Dict:
    """

    :param vals: the values, as a list or a numpy array
    :type vals: list or numpy array.
    :return: the statistics values.
    :rtype: a dictionary.
    """

    array = np.asarray(vals)

    min = np.nanmin(array)
    max = np.nanmax(array)
    mean = np.nanmean(array)
    var = np.nanvar(array)
    std = np.nanstd(array)

    stats = dict(min=min,
                 max=max,
                 mean=mean,
                 var=var,
                 std=std)

    return stats