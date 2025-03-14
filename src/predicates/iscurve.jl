# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    iscurve(g)

Tells whether or not the geometry `g` is a curve.
"""
iscurve(g::Geometry) = paramdim(g) == 1
iscurve(d::Domain) = paramdim(d) == 1
