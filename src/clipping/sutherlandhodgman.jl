# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SutherlandHodgman()

The Sutherland-Hodgman algorithm for clipping polygons.

## References

* Sutherland, I.E. & Hodgman, G.W. 1974. [Reentrant Polygon Clipping]
  (https://dl.acm.org/doi/pdf/10.1145/360767.360802)
"""
struct SutherlandHodgman <: ClippingMethod end

function clip(poly::Polygon, other::Polygon, ::SutherlandHodgman)
  v = vertices(poly)

  for s₁ in segments(other)
    newv = []
    n = length(v)

    for i in 1:n
      p₁, p₂ = v[i], v[i%n + 1]
      s₂ = Segment(p₁, p₂)

      # assuming convex clockwise other
      isp₁inside = (sideof(p₁, s₁) != :LEFT)
      isp₂inside = (sideof(p₂, s₁) != :LEFT)

      if isp₁inside && isp₂inside
        push!(newv, p₁)
      elseif isp₁inside && !isp₂inside
        push!(newv, p₁)
        push!(newv, s₁ ∩ s₂)
      elseif !isp₁inside && isp₂inside
        push!(newv, s₁ ∩ s₂)
      end
    end
    v = newv
  end
  poly = PolyArea(Ring(v...))
end
