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

"""
    intersecttype(s1, s2)

Compute the intersection type of two line segments `s1` and `s2` in 3D.
See https://stackoverflow.com/questions/55220355/how-to-detect-whether-two-segmentin-3d-spaceintersect.
"""
function intersecttype(s1::Segment{3,T}, s2::Segment{3,T}) where {T}
  # get coordinates
  p₁, p₂ = coordinates.(vertices(s1))
  q₁, q₂ = coordinates.(vertices(s2))

  # create matrices and vectors to allow for evaluation of intersection
  A = [(p₂ - p₁) (q₁ - q₂)]
  b = q₁ - p₁

  λ = A \ b

  # calculate the rank of the augmented matrix
  rₐ = rank([A b])
  # calculate the rank of the rectangular matrix
  r = rank(A)

  # use matrix rank to determine basic intersection properties
  # segments are co-planar (but not co-linear)
  if (rₐ == 2) && (r == 2)
    # if either element is approximately 0 or 1, set them as so to prevent any domain errors
    λ₁ = isapprox(λ[1], zero(T), atol=atol(T)) ? zero(T) : (isapprox(λ[1], one(T), atol=atol(T)) ? one(T) : λ[1])
    λ₂ = isapprox(λ[2], zero(T), atol=atol(T)) ? zero(T) : (isapprox(λ[2], one(T), atol=atol(T)) ? one(T) : λ[2])

    # if λs are outside of the interval [0, 1], they do not intersect
    if any((λ₁, λ₂) .< zero(T)) || any((λ₁, λ₂) .> one(T))
      NoIntersection()
    # if both λs are either 0 or 1, they are CornerTouchingSegments
    elseif ((λ₁ ≈ zero(T)) || (λ₁ ≈ one(T))) && ((λ₂ ≈ zero(T)) || (λ₂ ≈ one(T)))
      CornerTouchingSegments(s1(λ₁))
    # if either λ is 0 or 1 then they are MidTouchingSegments
    elseif ((λ₁ ≈ zero(T)) || (λ₁ ≈ one(T))) ⊻ ((λ₂ ≈ zero(T)) || (λ₂ ≈ one(T)))
      MidTouchingSegments(s1(λ₁))
    # otherwise they are simple CrossingSegments
    else
      CrossingSegments(s1(λ₁))
    end
  # segments are co-linear
  elseif (rₐ == 1) && (r == 1)
    intersectcolinear(s1, s2)
  else
    NoIntersection()
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

"""
    intersecttype(s, t)

Calculate the intersection of a segment and triangle in 3D.
Algorithm from Jiménez et al 2009 (A robust segment/triangle intersection algorithm for interference tests. Efficiency study). 
"""
function intersecttype(s::Segment{3,T}, t::Triangle{3,T}) where {T}
  vₛ = vertices(s)
  vₜ = vertices(t)
  
  A = vₛ[1] - vₜ[3]
  B = vₜ[1] - vₜ[3]
  C = vₜ[2] - vₜ[3]

  W₁ = B × C
  wᵥ = A ⋅ W₁

  D = vₛ[2] - vₜ[3]
  sᵥ = D ⋅ W₁

  if wᵥ > atol(T)
    # Rejection 2
    if sᵥ > atol(T)
      return NoIntersection()
    end

    W₂ = A × D
    tᵥ = W₂ ⋅ C
    
    # Rejection 3
    if tᵥ < -atol(T)
      return NoIntersection()
    end

    uᵥ = -(W₂ ⋅ B)

    # Rejection 4
    if uᵥ < -atol(T)
      return NoIntersection()
    end

    # Rejection 5
    if wᵥ < (sᵥ + tᵥ + uᵥ)
      return NoIntersection()
    end
  elseif wᵥ < -atol(T)
    # Rejection 2
    if sᵥ < -atol(T)
      return NoIntersection()
    end

    W₂ = A × D
    tᵥ = W₂ ⋅ C
    
    # Rejection 3
    if tᵥ > atol(T)
      return NoIntersection()
    end

    uᵥ = -(W₂ ⋅ B)

    # Rejection 4
    if uᵥ > atol(T)
      return NoIntersection()
    end

    # Rejection 5
    if wᵥ > (sᵥ + tᵥ + uᵥ)
      return NoIntersection()
    end
  else
    if sᵥ > atol(T)
      W₂ = D × A
      tᵥ = W₂ ⋅ C

      # Rejection 3
      if tᵥ < -atol(T)
        return NoIntersection()
      end

      uᵥ = -(W₂ ⋅ B)

      # Rejection 4
      if uᵥ < -atol(T)
        return NoIntersection()
      end
      # Rejection 5
      if -sᵥ < (tᵥ + uᵥ)
        return NoIntersection()
      end
    elseif sᵥ < -atol(T)
      W₂ = D × A
      tᵥ = W₂ ⋅ C

      # Rejection 3
      if tᵥ > atol(T)
        return NoIntersection()
      end

      uᵥ = -(W₂ ⋅ B)

      # Rejection 4
      if uᵥ > atol(T)
        return NoIntersection()
      end

      # Rejection 5
      if -sᵥ > (tᵥ + uᵥ)
        return NoIntersection()
      end
    else
      # Coplanar segment (Rejection 1, but we want to catch this)
      return intersectcoplanar(s, t)
    end
  end      

  λ = wᵥ / (wᵥ - sᵥ)

  # if λ is approximately 0 or 1, set as so to prevent any domain errors
  λ = isapprox(λ, zero(T), atol=atol(T)) ? zero(T) : (isapprox(λ, one(T), atol=atol(T)) ? one(T) : λ)

  # if λ is less than zero or greater than one, the triangle and segment don't intersect
  if (λ < zero(T)) || (λ > one(T))
    NoIntersection()
  else
    IntersectingSegmentTri(s(λ))
  end
end

"""
    intersectcoplanar(s, t)

Return the intersection of a Segment `s` and Triangle `t` given they are coplanar.

There are eight main cases to consider regarding `p₁`, `p₂`, the vertices making up the Segment `s` :

* `p₁`/`p₂` in `t`
  + `p₂`/`p₁` in `t`
    - 1. Returns a Segment
  + `p₂`/`p₁` on edge of `t`
    - 2. Returns a Segment 
  + `p₂`/`p₁` outside of `t`
    - 3. Returns a Segment
* `p₁`/`p₂` on edge of `t`
  + `p₂`/`p₁` on edge of `t`
    - 4. Returns a Segment 
  + `p₂`/`p₁` outside of `t`
    - 5. Returns a Point
* `p₁`/`p₂` outside of `t`
  + `p₂`/`p₁` outside of `t`
    - If the segment intersects an edge
      * 6. Returns a Segment
    - If the segment intersects a vertex
      * 7. Returns a Point
    - Otherwise
      * 8. Returns NoIntersection
"""
function intersectcoplanar(s::Segment{Dim,T}, t::Triangle{Dim,T}) where {Dim,T}
  vₛ = vertices(s)
  p₁, p₂ = vₛ

  # Determine segment intersections with triangle edges
  i = map(sₜ -> intersecttype(sₜ, s), segments(chains(t)[1]))

  # Get actual intersections and number of
  iᵥ = filter(i′ -> !isa(i′, NoIntersection), i)
  nᵢ = length(iᵥ)

  # Check for the special case if any intersection is an OverlappingSegments
  iₒ = filter(i′ -> isa(i′, OverlappingSegments), i)
  nₒ = length(iₒ)  
  
  # If both points of the segment are in the triangle, return the original Segment
  if (p₁ ∈ t) & (p₂ ∈ t)
    # Case 1, 2, 4
    OverlappingSegmentTri(s)
  elseif p₁ ∈ t
    # Segment is colinear with an edge
    if nₒ > 0
      # Case 5
      OverlappingSegmentTri(get(iₒ[1]))
    elseif nᵢ == 2
      # Get the edge intersections
      e₁ = get(iᵥ[1])
      e₂ = get(iᵥ[2])

      if e₁ ≈ e₂
        # Segment meets at a vertex of the triangle
        # Case 6
        IntersectingSegmentTri(e₁)
      else
        # Case 5
        OverlappingSegmentTri(Segment(e₁, e₂))
      end
    elseif nᵢ == 1
      # Get the edge intersection
      eᵢ = get(iᵥ[1])
      
      # If eᵢ and p₁ are the same, then the intersection is a point on an edge
      if eᵢ ≈ p₁
        # Case 6
        IntersectingSegmentTri(eᵢ)
      # Otherwise p₁ is "properly" inside the triangle and the intersection is a segment
      else 
        # Case 3
        OverlappingSegmentTri(Segment(p₁, eᵢ))
      end
    end
  elseif p₂ ∈ t
    # Segment is colinear with an edge
    if nₒ > 0
      # Case 5
      OverlappingSegmentTri(get(iₒ[1]))
    elseif nᵢ == 2
      # Get the edge intersections
      e₁ = get(iᵥ[1])
      e₂ = get(iᵥ[2])

      if e₁ ≈ e₂
        # Segment meets at a vertex of the triangle
        # Case 6
        IntersectingSegmentTri(e₁)      
      else
        # Case 5
        OverlappingSegmentTri(Segment(e₁, e₂))
      end
    elseif nᵢ == 1
      # Get the edge intersection
      eᵢ = get(iᵥ[1])
      
      # If eᵢ and p₂ are the same, then the intersection is a point on an edge
      if eᵢ ≈ p₂
        # Case 6
        IntersectingSegmentTri(eᵢ)
      # Otherwise p₂ is "properly" inside the triangle and the intersection is a segment
      else 
        # Case 3
        OverlappingSegmentTri(Segment(p₂, eᵢ))
      end
    end
  else
    # nᵢ can be either zero or two.
    if nᵢ == 0
      # Case 9
      NoIntersection()
    elseif nᵢ == 2
      # Get the edge intersections
      e₁ = get(iᵥ[1])
      e₂ = get(iᵥ[2])

      if e₁ ≈ e₂
        # Segment meets at a vertex of the triangle
        # Case 8
        IntersectingSegmentTri(e₁)
      else
        # Case 7
        OverlappingSegmentTri(Segment(e₁, e₂))
      end
    else
      # Shouldn't be possible, this serves as a "default"
      NoIntersection()
    end
  end
end
