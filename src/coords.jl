# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    coords(point)

Return the coordinates of the `point`.
"""
coords(A::Point) = A.coords

"""
    coords(vec)

Return the coordinates of the `vec`.
"""
coords(vec::StaticVector) = Cartesian(Tuple(vec))
