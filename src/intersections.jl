# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    g1 ∩ g2

Return the intersection of two geometries `g1` and `g2`.
"""
Base.intersect(g1::Geometry, g2::Geometry) = get(intersecttype(g1, g2))

"""
    intersecttype([f], g1, g2)

Compute the intersection type of two geometries `g1` and `g2`
and apply function `f` to it. Default function is [`identity`](@ref).

## Examples

```julia
intersecttype(g1, g2) do I
  if I isa CrossingLines
    # do something
  else
    # do nothing
  end
end
```

### Notes

When a custom function `f` is used that reduces the number of
return types, Julia is able to optimize the branches of the code
and generate specialized code. This is not the case when
`f === identity`.
"""
intersecttype(f, g1, g2) = intersecttype(f, g2, g1)
intersecttype(g1, g2)    = intersecttype(identity, g1, g2)

"""
    Intersection

An intersection type (e.g. crossing line segments, overlapping boxes).
"""
abstract type Intersection end

Base.get(I::Intersection) = I.value

# ------------------------
# LINE-LINE INTERSECTIONS
# ------------------------

struct CrossingLines{P<:Point} <: Intersection
  value::P
end

struct OverlappingLines{L<:Line} <: Intersection
  value::L
end

# ------------------------------
# SEGMENT-SEGMENT INTERSECTIONS
# ------------------------------

struct CrossingSegments{P<:Point} <: Intersection
  value::P
end

struct MidTouchingSegments{P<:Point} <: Intersection
  value::P
end

struct CornerTouchingSegments{P<:Point} <: Intersection
  value::P
end

struct OverlappingSegments{S<:Segment} <: Intersection
  value::S
end

# ----------------------
# BOX-BOX INTERSECTIONS
# ----------------------

struct OverlappingBoxes{B<:Box} <: Intersection
  value::B
end

struct FaceTouchingBoxes{B<:Box} <: Intersection
  value::B
end

struct CornerTouchingBoxes{P<:Point} <: Intersection
  value::P
end

# ----------------------
# RAY-BOX INTERSECTIONS
# ----------------------

struct CrossingRayBox{S<:Segment} <: Intersection
  value::S
end

struct TouchingRayBox{P<:Point} <: Intersection
  value::P
end

# -------------------------------
# SEGMENT-TRIANGLE INTERSECTIONS
# -------------------------------

struct IntersectingSegmentTriangle{P<:Point} <: Intersection
  value::P
end

# ---------------------------
# RAY-TRIANGLE INTERSECTIONS
# ---------------------------

struct IntersectingRayTriangle{P<:Point} <: Intersection
  value::P
end

# ----------------------------
# SEGMENT-PLANE INTERSECTIONS
# ----------------------------

struct CrossingSegmentPlane{P<:Point} <: Intersection
  value::P
end

struct TouchingSegmentPlane{P<:Point} <: Intersection
  value::P
end

struct OverlappingSegmentPlane{S<:Segment} <: Intersection
  value::S
end

# ------------
# CORNER CASE
# ------------

struct NoIntersection <: Intersection end

Base.get(::NoIntersection) = nothing

# ----------------
# IMPLEMENTATIONS
# ----------------

include("intersections/lines.jl")
include("intersections/segments.jl")
include("intersections/planes.jl")
include("intersections/triangles.jl")
include("intersections/boxes.jl")

"""
    hasintersect(g1, g2)

Return `true` if geometries `g1` and `g2` intersect and `false` otherwise.

The algorithm works with any geometry that has a well-defined [`supportfun`](@ref).

## References

* Gilbert, E., Johnson, D., Keerthi, S. 1988. [A fast
  Procedure for Computing the Distance Between Complex
  Objects in Three-Dimensional Space]
  (https://ieeexplore.ieee.org/document/2083)
"""
function hasintersect(g1::Geometry{Dim,T}, g2::Geometry{Dim,T}) where {Dim,T}
  # handle non-convex geometries
  if !isconvex(g1)
    d1 = discretize(g1, FIST())
    return hasintersect(d1, g2)
  elseif !isconvex(g2)
    d2 = discretize(g2, FIST())
    return hasintersect(g1, d2)
  end

  # initial direction
  c1, c2 = centroid(g1), centroid(g2)
  d = c1 == c2 ? rand(Vec{Dim,T}) : c2 - c1

  # first point in Minkowski difference
  P = minkowskipoint(g1, g2, d)

  # origin of coordinate system
  O = minkowskiorigin(Dim, T)

  # initialize simplex vertices
  points = [P]

  # move towards the origin
  d = O - P
  while true
    P = minkowskipoint(g1, g2, d)
    if (P - O) ⋅ d < zero(T)
      return false
    end
    push!(points, P)

    # line segment case
    if length(points) == 2
      B, A = points
      AB = B - A
      AO = O - A
      d  = AB × AO × AB
    else # simplex case
      C, B, A = points
      AB  = B - A
      AC  = C - A
      AO  = O - A
      ABᵀ = AC × AB × AB
      ACᵀ = AB × AC × AC
      if ABᵀ ⋅ AO > zero(T)
        popat!(points, 1) # pop C
        d = ABᵀ
      elseif ACᵀ ⋅ AO > zero(T)
        popat!(points, 2) # pop B
        d = ACᵀ
      else
        return true
      end
    end
  end
end

hasintersect(d1::Domain, g2::Geometry) =
  any(g1 -> hasintersect(g1, g2), d1)

hasintersect(g1::Geometry, d2::Domain) =
  hasintersect(d2, g1)

function hasintersect(d1::Domain, d2::Domain)
  for g1 in d1, g2 in d2
    hasintersect(g1, g2) && return true
  end
  return false
end

# support point in Minkowski difference
function minkowskipoint(g1::Geometry{Dim,T}, g2::Geometry{Dim,T}, d) where {Dim,T}
  n = Vec{Dim}(d[1:Dim])
  v = supportfun(g1, n) - supportfun(g2, -n)
  Point(ntuple(i-> i ≤ Dim ? v[i] : zero(T), max(Dim, 3)))
end

# origin of coordinate system
minkowskiorigin(Dim, T) = Point(ntuple(i->zero(T), max(Dim, 3)))
