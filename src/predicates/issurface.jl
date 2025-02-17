# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    issurface(g)

Tells whether or not the geometry `g` is a surface.
"""
issurface(g::Geometry) = paramdim(g) == 2
issurface(d::Domain) = paramdim(d) == 2
