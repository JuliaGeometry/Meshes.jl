# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    intersecttype(s, t)

Calculate the intersection of a segment and triangle in 3D.

## References

* Jiménez, J., Segura, R. and Feito, F. 2009. [A robust segment/triangle
  intersection algorithm for interference tests. Efficiency study]
  (https://www.sciencedirect.com/science/article/pii/S0925772109001448?via%3Dihub)
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
    # rejection 2
    if sᵥ > atol(T)
      return NoIntersection()
    end

    W₂ = A × D
    tᵥ = W₂ ⋅ C
    
    # rejection 3
    if tᵥ < -atol(T)
      return NoIntersection()
    end

    uᵥ = -(W₂ ⋅ B)

    # rejection 4
    if uᵥ < -atol(T)
      return NoIntersection()
    end

    # rejection 5
    if wᵥ < (sᵥ + tᵥ + uᵥ)
      return NoIntersection()
    end
  elseif wᵥ < -atol(T)
    # rejection 2
    if sᵥ < -atol(T)
      return NoIntersection()
    end

    W₂ = A × D
    tᵥ = W₂ ⋅ C
    
    # rejection 3
    if tᵥ > atol(T)
      return NoIntersection()
    end

    uᵥ = -(W₂ ⋅ B)

    # rejection 4
    if uᵥ > atol(T)
      return NoIntersection()
    end

    # rejection 5
    if wᵥ > (sᵥ + tᵥ + uᵥ)
      return NoIntersection()
    end
  else
    if sᵥ > atol(T)
      W₂ = D × A
      tᵥ = W₂ ⋅ C

      # rejection 3
      if tᵥ < -atol(T)
        return NoIntersection()
      end

      uᵥ = -(W₂ ⋅ B)

      # rejection 4
      if uᵥ < -atol(T)
        return NoIntersection()
      end
      # rejection 5
      if -sᵥ < (tᵥ + uᵥ)
        return NoIntersection()
      end
    elseif sᵥ < -atol(T)
      W₂ = D × A
      tᵥ = W₂ ⋅ C

      # rejection 3
      if tᵥ > atol(T)
        return NoIntersection()
      end

      uᵥ = -(W₂ ⋅ B)

      # rejection 4
      if uᵥ > atol(T)
        return NoIntersection()
      end

      # rejection 5
      if -sᵥ > (tᵥ + uᵥ)
        return NoIntersection()
      end
    else
      # rejection 1, coplanar segment
      return NoIntersection()
    end
  end      

  λ = wᵥ / (wᵥ - sᵥ)

  # if λ is approximately 0 or 1, set as so to prevent any domain errors
  λ = isapprox(λ, zero(T), atol=atol(T)) ? zero(T) : (isapprox(λ, one(T), atol=atol(T)) ? one(T) : λ)
  
  IntersectingSegmentTriangle(s(λ))
end