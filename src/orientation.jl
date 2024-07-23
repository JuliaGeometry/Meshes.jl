# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    OrientationType

The different types of orientation of a ring.
Possible values are `CW` and `CCW`.
"""
@enum OrientationType begin
  CW
  CCW
end

"""
    orientation(geom)

Returns the orientation of the geometry `geom` as
either counter-clockwise (CCW) or clockwise (CW).
"""
function orientation end

function orientation(p::Polygon)
  o = [orientation(ring) for ring in rings(p)]
  hasholes(p) ? o : first(o)
end

orientation(r::Ring) = _orientation(r, Val(embeddim(r)))

_orientation(r::Ring, ::Val{3}) = _orientation(proj2D(r), Val(2))

function _orientation(r::Ring, ::Val{2})
  ℒ = lentype(r)
  v = vertices(r)
  n = nvertices(r)
  A(i) = signarea(v[1], v[i], v[i + 1])
  Σ = sum(A, 2:(n - 1), init=zero(ℒ)^2)
  Σ ≥ zero(Σ) ? CCW : CW
end
