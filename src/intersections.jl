# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    g1 ∩ g2

Return the intersection of two geometries `g1` and `g2` as a new geometry.
"""
Base.intersect(g1::Geometry, g2::Geometry) = get(intersection(g1, g2))

"""
    intersection([f], g1, g2)

Compute the intersection of two geometries `g1` and `g2`
and apply function `f` to it. Default function is `identity`.

## Examples

```julia
intersection(g1, g2) do I
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
intersection(f, g1, g2) = intersection(f, g2, g1)
intersection(g1, g2) = intersection(identity, g1, g2)

"""
    IntersectionType

The different types of intersection that may occur between geometries.
Type `IntersectionType` in a Julia session to see the full list.
"""
@enum IntersectionType begin
  # no intersection
  NoIntersection

  # line-line intersection
  CrossingLines
  OverlappingLines

  # box-box intersection
  OverlappingBoxes
  FaceTouchingBoxes
  CornerTouchingBoxes

  # segment-segment intersection
  CrossingSegments
  MidTouchingSegments
  CornerTouchingSegments
  OverlappingSegments

  # ray-ray intersection
  CrossingRays
  MidTouchingRays
  CornerTouchingRays
  OverlappingAgreeingRays
  OverlappingOpposingRays

  # ray-segment intersection
  CrossingRaySegment
  MidTouchingRaySegment
  CornerTouchingRaySegment
  OverlappingRaySegment

  # ray-line intersection
  CrossingRayLine
  TouchingRayLine
  OverlappingRayLine

  # line-segment intersection
  CrossingLineSegment
  TouchingLineSegment
  OverlappingLineSegment

  # ray-box intersection
  CrossingRayBox
  TouchingRayBox

  # segment-triangle intersection
  IntersectingSegmentTriangle

  # ray-triangle intersection
  CrossingRayTriangle
  TouchingRayTriangle
  CornerTouchingRayTriangle
  CornerCrossingRayTriangle
  EdgeCrossingRayTriangle
  EdgeTouchingRayTriangle

  # line-plane intersection
  CrossingLinePlane
  OverlappingLinePlane

  # ray-plane intersection
  CrossingRayPlane
  TouchingRayPlane
  OverlappingRayPlane

  # segment-plane intersection
  CrossingSegmentPlane
  TouchingSegmentPlane
  OverlappingSegmentPlane
end

"""
    Intersection{G}

An intersection between geometries holding a geometry of type `G`.
"""
struct Intersection{GeometryType}
  type::IntersectionType
  geom::GeometryType
end

Intersection(type, geom) = Intersection{typeof(geom)}(type, geom)

"""
    type(intersection)

Return the type of intersection computed between geometries.
"""
type(I::Intersection) = I.type

"""
    get(intersection)

Return the underlying geometry stored in the intersection object.
"""
Base.get(I::Intersection) = I.geom

# helper macro for developers in case we decide to
# change the internal representation of Intersection
macro IT(type, geom, func)
  type = esc(type)
  geom = esc(geom)
  func = esc(func)
  :(Intersection($type, $geom) |> $func)
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("intersections/lines.jl")
include("intersections/boxes.jl")
include("intersections/segments.jl")
include("intersections/rays.jl")
include("intersections/raysegment.jl")
include("intersections/rayline.jl")
include("intersections/linesegment.jl")
include("intersections/raybox.jl")
include("intersections/lineraysegmentplane.jl")
include("intersections/segmenttriangle.jl")
include("intersections/raytriangle.jl")
include("intersections/geompolygon.jl")

# ------------------------
# MISSING IMPLEMENTATIONS
# ------------------------

intersection(f, ::Ball, ::Polygon) = throw(ErrorException("not implemented"))

# -----------------
# TRUE/FALSE CHECK
# -----------------

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
    d1 = simplexify(g1)
    return hasintersect(d1, g2)
  elseif !isconvex(g2)
    d2 = simplexify(g2)
    return hasintersect(g1, d2)
  end

  # initial direction
  c1, c2 = centroid(g1), centroid(g2)
  d = c1 ≈ c2 ? rand(Vec{Dim,T}) : c2 - c1

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
      d = AB × AO × AB
    else # simplex case
      C, B, A = points
      AB = B - A
      AC = C - A
      AO = O - A
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

hasintersect(p::Point, g::Geometry) = p ∈ g

hasintersect(g::Geometry, p::Point) = hasintersect(p, g)

hasintersect(p1::Point, p2::Point) = p1 == p2

hasintersect(d::Domain, g::Geometry) = any(gᵢ -> hasintersect(gᵢ, g), d)

hasintersect(g::Geometry, d::Domain) = hasintersect(d, g)

function hasintersect(d1::Domain, d2::Domain)
  for g1 in d1, g2 in d2
    hasintersect(g1, g2) && return true
  end
  return false
end

hasintersect(g) = Base.Fix2(hasintersect, g)

# support point in Minkowski difference
function minkowskipoint(g1::Geometry{Dim,T}, g2::Geometry{Dim,T}, d) where {Dim,T}
  n = Vec{Dim,T}(d[1:Dim])
  v = supportfun(g1, n) - supportfun(g2, -n)
  Point(ntuple(i -> i ≤ Dim ? v[i] : zero(T), max(Dim, 3)))
end

# origin of coordinate system
minkowskiorigin(Dim, T) = Point(ntuple(i -> zero(T), max(Dim, 3)))
