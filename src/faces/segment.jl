# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# ORDERING CONVENTION
#       v
#       ^
#       |
#       |
# 0-----+-----1 --> u

"""
    Segment(p1, p2)

A line segment with points `p1`, `p2`.
"""
struct Segment{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Face{Dim,T}
  vertices::V
end

paramdim(::Type{<:Segment}) = 1

facets(s::Segment) = (v for v in s.vertices)

"""
    s1 ∩ s2

Compute the intersection of two line segments `s1` and `s2`.
"""
Base.intersect(s1::Segment, s2::Segment) = get(intersecttype(s1, s2))

abstract type Intersection end

Base.get(I::Intersection) = I.value

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

struct NonIntersectingSegments <: Intersection end

Base.get(::NonIntersectingSegments) = nothing

"""
    intersecttype(s1, s2)

Compute the intersection type of two line segments `s1` and `s2`.

The intersection type can be one of five types according to
Balbes, R. and Siegel, J. 1990:

1. intersect at point interior to both segments
2. intersect at end point of one and an interior point of the other
3. intersect at end point of both segments
4. overlap at more than one point
5. do not overlap nor intersect

## References

* Balbes, R. and Siegel, J. 1990. [A robust method for calculating
  the simplicity and orientation of planar polygons]
  (https://www.sciencedirect.com/science/article/abs/pii/0167839691900198)
"""
function intersecttype(s1::Segment{2,T}, s2::Segment{2,T}) where {T}
  x̄ = center(s1)
  ȳ = center(s2)
  x1, x2 = vertices(s1)
  y1, y2 = vertices(s2)
  Q = [x1, y1, x2, y2, x1]

  # winding number of x̄
  windingx = zero(T)
  determinatex = true
  for j in 1:length(Q)-1
    α = ∠(Q[j], x̄, Q[j+1])
    if α ≈ 0 || α ≈ π || α ≈ -π
      determinatex = false
      break
    end
    windingx += α
  end

  # winding number of ȳ
  windingy = zero(T)
  determinatey = true
  for j in 1:length(Q)-1
    β = ∠(Q[j], ȳ, Q[j+1])
    if β ≈ 0 || β ≈ π || β ≈ -π
      determinatey = false
      break
    end
    windingy += β
  end

  # winding number ≈ 0?
  wxzero = isapprox(windingx, 0, atol=atol(T))
  wyzero = isapprox(windingy, 0, atol=atol(T))

  # there are three cases to consider
  # Balbes, R. and Siegel, J. 1990.
  if determinatex && determinatey # CASE (I)
    if !wxzero && !wyzero
      # configuration (1)
      CrossingSegments(intersectpoint(s1, s2))
    else
      # configuration (5)
      NonIntersectingSegments()
    end
  elseif determinatex || determinatey # CASE (II)
    if !(determinatex ? wxzero : wyzero)
      if x1 == y1 || x1 == y2 || x2 == y1 || x2 == y2
        # configuration (3)
        CornerTouchingSegments(intersectpoint(s1, s2))
      else
        # configuration (2)
        MidTouchingSegments(intersectpoint(s1, s2))
      end
    else
      # configuration (5)
      NonIntersectingSegments()
    end
  elseif !determinatex && !determinatey # CASE (III)
    if !isapprox((x2 - x1) × (y2 - y1), 0, atol=atol(T))
      # configuration (3)
      CornerTouchingSegments(intersectpoint(s1, s2))
    else
      if x1 == y1 || x1 == y2 || x2 == y1 || x2 == y2
        # configuration (3)
        CornerTouchingSegments(intersectpoint(s1, s2))
      elseif x2 ≠ y1
        # configuration (4)
        OverlappingSegments(intersectsegment(s1, s2))
      else
        # configuration (5)
        NonIntersectingSegments()
      end
    end
  end
end

"""
    intersectpoint(s1, s2)

Compute the intersection of two line segments `s1` and `s2`
assuming that it is a point.
"""
function intersectpoint(s1::Segment{Dim,T}, s2::Segment{Dim,T}) where {Dim,T}
  x1, x2 = vertices(s1)
  y1, y2 = vertices(s2)

  # solve Ax = b
  A = (y2 - y1) - (x2 - x1)
  b = (x1 - y1)
  s = zero(T)
  for i in 1:length(A)
    if !iszero(A[i])
      s = b[i] / A[i]
      break
    end
  end

  x1 + s*(x2 - x1)
end

"""
    intersectsegment(s1, s2)

Compute the intersection of two line segments `s1` and `s2`
assuming that it is a segment.
"""
function intersectsegment(s1::Segment{Dim,T}, s2::Segment{Dim,T}) where {Dim,T}
  x1, x2 = vertices(s1)
  y1, y2 = vertices(s2)

  # solve Ax = b
  A = x2 - x1
  b = y1 - x1
  s = zero(T)
  for i in 1:length(A)
    if !iszero(A[i])
      s = b[i] / A[i]
      break
    end
  end

  if s > 0 # s1 is to the left of s2
    Segment(y1, x2)
  else # s1 is to the right of s2
    Segment(x1, y2)
  end
end
