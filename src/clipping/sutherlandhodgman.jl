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
    l₁ = Line(vertices(s₁)...)
    n = length(v)
    newv = []
    
    for i in 1:n
      p₁, p₂ = v[i], v[i%n + 1]
      l₂ = Line(p₁, p₂)

      # assuming convex clockwise other
      isp₁inside = (sideof(p₁, l₁) != :LEFT)
      isp₂inside = (sideof(p₂, l₁) != :LEFT)

      if isp₁inside && isp₂inside
        push!(newv, p₁)
      elseif isp₁inside && !isp₂inside
        push!(newv, p₁)
        push!(newv, l₁ ∩ l₂)
      elseif !isp₁inside && isp₂inside
        push!(newv, l₁ ∩ l₂)
      end
    end
    v = newv
  end
  poly = PolyArea(Ring(v...))
end
