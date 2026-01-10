# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
  MartinezRuedaClipping()

The Martinez-Rueda algorithm for clipping arbitrary polygons.

## References

* Martínez, F., Rueda, A.J., Feito, F.R. 2009. [A new algorithm for computing Boolean operations on
  polygons](https://doi.org/10.1016/j.cag.2009.03.003)
"""
struct MartinezRuedaClipping <: ClippingMethod end

function clip(poly::Polygon, other::Geometry, method::MartinezRuedaClipping)
  polygonbooleanop(poly, other, intersect)
end

function clip(ring::Ring, other::Ring, ::MartinezRuedaClipping)
  polygonbooleanop(ring, other, intersect)
end
