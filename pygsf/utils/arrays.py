import numbers
from array import array
from typing import List

from pygsf.utils.types import check_type


class ArrayList:
    """
    An identified list of arrays.
    """

    def __init__(self,
                 rec_id: numbers.Integral,
                 arrays: List[array]
                 ):

        check_type(rec_id, "Input id", numbers.Integral)
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

