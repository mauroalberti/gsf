import numbers
from array import array
from typing import Union, List, Tuple

from numpy import linspace

from pygsf.utils.types import check_type


class ArrayList:
    """
    An identified list of arrays.
    """

    def __init__(self,
                 rec_id: Union[str, numbers.Integral],
                 arrays: List[array]
                 ):

        check_type(rec_id, "Input id", (str, numbers.Integral))
        check_type(arrays, "List of arrays", list)
        for part in arrays:
            check_type(part, "Array", array)

        self._id = rec_id
        self._arrays = arrays

    @property
    def id(self) -> numbers.Integral:
        """
        Return id.

        :return: the id
        :rtype: numbers.Integral
        """

        return self._id

    @property
    def arrays(self) -> List[array]:
        """
        Returns the object arrays.

        :return: the instance arrays
        :rtype: List[array]
        """

        return self._arrays


def to_float(
        curr_iterable
) -> Tuple[numbers.Real]:

    return (float(item) for item in curr_iterable)


def almost_zero(val):

    tolerance = 1e-10
    if abs(val) > tolerance:
        return False
    else:
        return True


def formula_to_grid(array_range, array_size, formula):

    a_min, a_max, b_max, b_min = array_range  # note: b range reversed for conventional j order in arrays
    array_rows, array_cols = array_size

    a_array = linspace(a_min, a_max, num=array_cols)
    b_array = linspace(b_max, b_min, num=array_rows)  # note: reversed for conventional j order in arrays

    try:
        a_list, b_list = [a for a in a_array for b in b_array], [b for a in a_array for b in b_array]
    except:
        raise Exception("Error in a-b values")

    try:
        z_list = [eval(formula) for a in a_array for b in b_array]
    except:
        raise Exception("Error in applying formula to a and b array values")

    return a_list, b_list, z_list

