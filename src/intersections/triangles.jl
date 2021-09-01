# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

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
    IntersectingSegmentTriangle(s(λ))
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
    OverlappingSegmentTriangle(s)
  elseif p₁ ∈ t
    # Segment is colinear with an edge
    if nₒ > 0
      # Case 5
      OverlappingSegmentTriangle(get(iₒ[1]))
    elseif nᵢ == 2
      # Get the edge intersections
      e₁ = get(iᵥ[1])
      e₂ = get(iᵥ[2])

      if e₁ ≈ e₂
        # Segment meets at a vertex of the triangle
        # Case 6
        IntersectingSegmentTriangle(e₁)
      else
        # Case 5
        OverlappingSegmentTriangle(Segment(e₁, e₂))
      end
    elseif nᵢ == 1
      # Get the edge intersection
      eᵢ = get(iᵥ[1])
      
      # If eᵢ and p₁ are the same, then the intersection is a point on an edge
      if eᵢ ≈ p₁
        # Case 6
        IntersectingSegmentTriangle(eᵢ)
      # Otherwise p₁ is "properly" inside the triangle and the intersection is a segment
      else 
        # Case 3
        OverlappingSegmentTriangle(Segment(p₁, eᵢ))
      end
    end
  elseif p₂ ∈ t
    # Segment is colinear with an edge
    if nₒ > 0
      # Case 5
      OverlappingSegmentTriangle(get(iₒ[1]))
    elseif nᵢ == 2
      # Get the edge intersections
      e₁ = get(iᵥ[1])
      e₂ = get(iᵥ[2])

      if e₁ ≈ e₂
        # Segment meets at a vertex of the triangle
        # Case 6
        IntersectingSegmentTriangle(e₁)      
      else
        # Case 5
        OverlappingSegmentTriangle(Segment(e₁, e₂))
      end
    elseif nᵢ == 1
      # Get the edge intersection
      eᵢ = get(iᵥ[1])
      
      # If eᵢ and p₂ are the same, then the intersection is a point on an edge
      if eᵢ ≈ p₂
        # Case 6
        IntersectingSegmentTriangle(eᵢ)
      # Otherwise p₂ is "properly" inside the triangle and the intersection is a segment
      else 
        # Case 3
        OverlappingSegmentTriangle(Segment(p₂, eᵢ))
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
        IntersectingSegmentTriangle(e₁)
      else
        # Case 7
        OverlappingSegmentTriangle(Segment(e₁, e₂))
      end
    else
      # Shouldn't be possible, this serves as a "default"
      NoIntersection()
    end
  end
end