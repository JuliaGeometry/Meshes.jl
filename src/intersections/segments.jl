# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

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
  x̄ = centroid(s1)
  ȳ = centroid(s2)
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
  wxzero = isapprox(windingx, zero(T), atol=atol(T))
  wyzero = isapprox(windingy, zero(T), atol=atol(T))

  # there are three cases to consider
  # Balbes, R. and Siegel, J. 1990.
  if determinatex && determinatey # CASE (I)
    if !wxzero && !wyzero
      # configuration (1)
      CrossingSegments(intersectpoint(s1, s2))
    else
      # configuration (5)
      NoIntersection()
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
      NoIntersection()
    end
  elseif !determinatex && !determinatey # CASE (III)
    if !isapprox((x2 - x1) × (y2 - y1), zero(T), atol=atol(T)^2)
      # configuration (3)
      CornerTouchingSegments(intersectpoint(s1, s2))
    else
      # configuration (3), (4) or (5)
      intersectcolinear(s1, s2)
    end
  end
end

# compute the intersection of two line segments assuming that it is a point
function intersectpoint(s1::Segment{2}, s2::Segment{2})
  x1, x2 = vertices(s1)
  y1, y2 = vertices(s2)
  intersectpoint(Line(x1, x2), Line(y1, y2))
end

# intersection of two line segments assuming that they are colinear
function intersectcolinear(s1::Segment{Dim,T}, s2::Segment{Dim,T}) where {Dim,T}
  m1, M1 = coordinates.(vertices(s1))
  m2, M2 = coordinates.(vertices(s2))

  # make sure that segment vertices are "ordered"
  m1, M1 = any(m1 .> M1) ? (M1, m1) : (m1, M1)
  m2, M2 = any(m2 .> M2) ? (M2, m2) : (m2, M2)

  # relevant vertices
  u = Point(max.(m1, m2))
  v = Point(min.(M1, M2))

  if isapprox(u, v, atol=atol(T))
    CornerTouchingSegments(u)
  elseif any(coordinates(u) .< coordinates(v))
    OverlappingSegments(Segment(u, v))
  else
    NoIntersection()
  end
end
