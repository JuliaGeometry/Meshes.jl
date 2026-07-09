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

orientation(r::Ring{𝔼{3}}) = orientation(proj2D(r))

function orientation(r::Ring{𝔼{2}})
  A = signedenclosedarea(r)
  A ≥ zero(A) ? CCW : CW
end

function orientation(r::Ring)
  ℒ = lentype(r)
  v = vertices(r)
  n = nvertices(r)
  A(i) = signarea(flat(v[1]), flat(v[i]), flat(v[i + 1]))
  Σ = sum(A, 2:(n - 1), init=zero(ℒ)^2)
  Σ ≥ zero(Σ) ? CCW : CW
end
