from .geometries import *


def try_derive_bestfitplane(
    points: Points
) -> Tuple[bool, Union[str, Plane]]:

    print(points.xs)
    print(points.ys)
    print(points.zs)

    npaXyz = points.asXyzArray()

    print(points.asXyzArray())

    xyz_mean = np.mean(npaXyz, axis=0)

    svd = xyzSvd(npaXyz - xyz_mean)

    if svd['result'] is None:
        return False, "Unable to calculate result"

    _, _, eigenvectors = svd['result']

    lowest_eigenvector = eigenvectors[-1, : ]  # Solution is last row

    normal = lowest_eigenvector[: 3 ] / np.linalg.norm(lowest_eigenvector[: 3 ])
    normal_vector = Vect(normal[0], normal[1], normal[2])
    normal_direct = Direct.fromVect(normal_vector)

    return True, normal_direct.normPlane()

