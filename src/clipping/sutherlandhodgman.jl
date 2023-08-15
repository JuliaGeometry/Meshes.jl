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

function clip(poly::Polygon{Dim,T}, other::Polygon, ::SutherlandHodgman) where {Dim,T}
  v = vertices(poly)

  for s₁ in segments(other)
    l₁ = Line(vertices(s₁)...)
    n = length(v)
    newv::Vector{Point{Dim,T}} = []
    
    for i in 1:n
      p₁, p₂ = v[i], v[mod1(i+1, n)]
      l₂ = Line(p₁, p₂)

      # assuming convex clockwise other
      isinside₁ = (sideof(p₁, l₁) != :LEFT)
      isinside₂ = (sideof(p₂, l₁) != :LEFT)

      if isinside₁ && isinside₂
        push!(newv, p₁)
      elseif isinside₁ && !isinside₂
        push!(newv, p₁)
        push!(newv, l₁ ∩ l₂)
      elseif !isinside₁ && isinside₂
        push!(newv, l₁ ∩ l₂)
      end
    end
    v = newv
  end
  poly = PolyArea(Ring(v))
end
