# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    issolid(g)

Tells whether or not the geometry `g` is a solid (i.e., 3D volume).
"""
issolid(g::Geometry) = paramdim(g) == 3
issolid(d::Domain) = paramdim(d) == 3
