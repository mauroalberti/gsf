

def trend(self) -> Optional[float]:
    """
    Trend of a vector
    (degrees, clockwise from North, range 0°-360°)

    :return: an optional float value, representing trend, in decimal degrees.

    Examples:
      >>> trend(1, 0, 0)
      90.0
      >>> trend(0, 1, 0)
      0.0
      >>> trend(1, 1, 0)
      45.0
      >>> trend(1, -1, 0)
      135.0
      >>> trend(0, -1, 0)
      180.0
      >>> trend(-1, -1, 0)
      225.0
      >>> trend(-1, 0, 0)
      270.0
      >>> trend(-1, 1, 0)
      315.0
      >>> trend(1, 1, 10)
      45.0
      >>> trend(0, 0, 0)
      None
     """

    return angle_north_clock(self.x. self.y)

def slope(self) -> Optional[float]:
    """
    Slope of a vector.
    Degrees, positive: downward-directed, negative: upward-dir., range -90°/90°

    :return: an optional float, representing the vector slope, in decimal degrees.

    Examples:
      >>> slope(Vect(1, 0, -1)
      45.0
      >>> slope(Vect(1, 0, 1)
      (-45.0)
      >>> slope(Vect(0, 1, 0)
      0.0
      >>> slope(Vect(0, 0, 1)
      (-90.0)
      >>> slope(Vect(0, 0, -1)
      90.0
      >>> slope(Vect(0, 0, 0)
      Nothing
     """
    h = v.len2D()
    zv = self.z
    sl = slope(h, abs(zv))

    if sl is None:
        return None
    else:
        if zv <= 0.0:
            return sl
        else:
            return -sl


