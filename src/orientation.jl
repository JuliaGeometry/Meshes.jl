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
    OrientationMethod

A method for finding the orientation of rings and polygons.
"""
abstract type OrientationMethod end

"""
    orientation(geom, [method])

Returns the orientation of the geometry `geom` as
either counter-clockwise (CCW) or clockwise (CW).

Optionally, specify the orientation `method`.

See also [`WindingOrientation`](@ref),
[`TriangleOrientation`](@ref).
"""
function orientation end

orientation(p::Polygon) = orientation(p, WindingOrientation())

orientation(r::Ring) = orientation(r, WindingOrientation())

function orientation(p::Polygon, method)
  o = [orientation(ring, method) for ring in rings(p)]
  hasholes(p) ? o : first(o)
end

orientation(r::Ring{3}, method) = orientation(proj2D(r), method)

"""
    WindingOrientation()

A method for finding the orientatino of rings and polygons
based on the winding number.

## References

* Balbes, R. and Siegel, J. 1990. [A robust method for calculating
  the simplicity and orientation of planar polygons]
  (https://www.sciencedirect.com/science/article/abs/pii/0167839691900198)
"""
struct WindingOrientation <: OrientationMethod end

function orientation(r::Ring{2}, ::WindingOrientation)
  # pick any segment
  p₁, p₂ = vertices(r)[1:2]
  p̄ = center(Segment(p₁, p₂))
  w̄ = winding(p̄, r)
  w = oftype(w̄, 2π) * w̄ - ∠(p₁, p̄, p₂)
  isapproxequal(w, oftype(w, π)) ? CCW : CW
end

"""
    TriangleOrientation()

A method for finding the orientation of rings and polygons
based on signed triangular areas.

## References

* Held, M. 1998. [FIST: Fast Industrial-Strength Triangulation of Polygons]
  (https://link.springer.com/article/10.1007/s00453-001-0028-4)
"""
struct TriangleOrientation <: OrientationMethod end

function orientation(r::Ring{2}, ::TriangleOrientation)
  v = vertices(r)
  Δ(i) = signarea(v[1], v[i], v[i + 1])
  a = mapreduce(Δ, +, 2:(length(v) - 1))
  a ≥ zero(a) ? CCW : CW
end
