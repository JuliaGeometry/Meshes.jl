# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ClippingMethod

A method for clipping geometries with other geometries.
"""
abstract type ClippingMethod end

"""
    clip(geometry, other, method)
"""
function clip(poly::Polygon, other::Geometry, method::ClippingMethod)
  c = [clip(ring, boundary(other), method) for ring in rings(poly)]
  r = filter(!isnothing, c)

  if isempty(r)
    nothing
  else
    PolyArea(r)
  end
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("clipping/sutherlandhodgman.jl")
