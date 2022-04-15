# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    intersection(segment, triangle)

Compute the intersection type of `segment and `triangle`.

## References

* Jiménez, J., Segura, R. and Feito, F. 2009. [A robust segment/triangle
  intersection algorithm for interference tests. Efficiency study]
  (https://www.sciencedirect.com/science/article/pii/S0925772109001448?via%3Dihub)
"""
function intersection(s::Segment{3,T}, t::Triangle{3,T}) where {T}
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
      return @IT NoIntersection nothing
    end

    W₂ = A × D
    tᵥ = W₂ ⋅ C

    # rejection 3
    if tᵥ < -atol(T)
      return @IT NoIntersection nothing
    end

    uᵥ = -(W₂ ⋅ B)

    # rejection 4
    if uᵥ < -atol(T)
      return @IT NoIntersection nothing
    end

    # rejection 5
    if wᵥ < (sᵥ + tᵥ + uᵥ)
      return @IT NoIntersection nothing
    end
  elseif wᵥ < -atol(T)
    # rejection 2
    if sᵥ < -atol(T)
      return @IT NoIntersection nothing
    end

    W₂ = A × D
    tᵥ = W₂ ⋅ C

    # rejection 3
    if tᵥ > atol(T)
      return @IT NoIntersection nothing
    end

    uᵥ = -(W₂ ⋅ B)

    # rejection 4
    if uᵥ > atol(T)
      return @IT NoIntersection nothing
    end

    # rejection 5
    if wᵥ > (sᵥ + tᵥ + uᵥ)
      return @IT NoIntersection nothing
    end
  else
    if sᵥ > atol(T)
      W₂ = D × A
      tᵥ = W₂ ⋅ C

      # rejection 3
      if tᵥ < -atol(T)
        return @IT NoIntersection nothing
      end

      uᵥ = -(W₂ ⋅ B)

      # rejection 4
      if uᵥ < -atol(T)
        return @IT NoIntersection nothing
      end
      # rejection 5
      if -sᵥ < (tᵥ + uᵥ)
        return @IT NoIntersection nothing
      end
    elseif sᵥ < -atol(T)
      W₂ = D × A
      tᵥ = W₂ ⋅ C

      # rejection 3
      if tᵥ > atol(T)
        return @IT NoIntersection nothing
      end

      uᵥ = -(W₂ ⋅ B)

      # rejection 4
      if uᵥ > atol(T)
        return @IT NoIntersection nothing
      end

      # rejection 5
      if -sᵥ > (tᵥ + uᵥ)
        return @IT NoIntersection nothing
      end
    else
      # rejection 1, coplanar segment
      return @IT NoIntersection nothing
    end
  end

  λ = clamp(wᵥ / (wᵥ - sᵥ), zero(T), one(T))

  return @IT IntersectingSegmentTriangle s(λ)
end
