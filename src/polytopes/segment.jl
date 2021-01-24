# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Segment(p1, p2)

A line segment with points `p1`, `p2`.
"""
struct Segment{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Polytope{1,Dim,T}
  vertices::V
end

measure(s::Segment) = norm(s.vertices[2] - s.vertices[1])

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
      if x1 ≈ y1 || x1 ≈ y2 || x2 ≈ y1 || x2 ≈ y2
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
      # configuration (3), (4) or (5)
      intersectcolinear(s1, s2)
    end
  end
end

# compute the intersection of two line segments assuming that it is a point
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

# intersection of two line segments assuming that they are colinear
function intersectcolinear(s1::Segment{Dim,T}, s2::Segment{Dim,T}) where {Dim,T}
  # make sure the first segment is larger than the second
  sa, sb = measure(s1) < measure(s2) ? (s2, s1) : (s1, s2)
  x1, x2 = vertices(sa)
  y1, y2 = vertices(sb)

  # given that the segments are colinear we can simply
  # operate on the scalar parameters of the vertices
  # along the line defined by the vector --x1--x2->
  Δ = x2 - x1
  i = findfirst(!iszero, Δ)
  t1 = (y1 - x1)[i] / Δ[i]
  t2 = (y2 - x1)[i] / Δ[i]

  # fix direction of second segment to match the first
  t1 > t2 && ((t1, t2) = (t2, t1))

  if t1 < zero(T) && t2 < zero(T)
    NonIntersectingSegments()
  elseif t1 < zero(T) && t2 ≈ zero(T)
    CornerTouchingSegments(x1)
  elseif t1 < zero(T) && t2 > zero(T)
    OverlappingSegments(Segment(x1, x1 + t2*Δ))
  elseif t1 ≈ zero(T) && t2 > zero(T)
    OverlappingSegments(Segment(x1, x1 + t2*Δ))
  elseif t1 > zero(T) && t2 < one(T)
    OverlappingSegments(Segment(x1 + t1*Δ, x1 + t2*Δ))
  elseif t1 < one(T) && t2 ≈ one(T)
    OverlappingSegments(Segment(x1 + t1*Δ, x2))
  elseif t1 < one(T) && t2 > one(T)
    OverlappingSegments(Segment(x1 + t1*Δ, x1 + t2*Δ))
  elseif t1 ≈ one(T) && t2 > one(T)
    CornerTouchingSegments(x2)
  elseif t1 > one(T) && t2 > one(T)
    NonIntersectingSegments()
  else
    @error "please report bug"
  end
end
