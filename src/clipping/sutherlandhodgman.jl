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

function clip(poly::Polygon{Dim,T}, other::Polygon{Dim,T}, ::SutherlandHodgman) where {Dim,T}
  vother = vertices(other)
  n = length(vother)
  v = vertices(poly)

  for i in 1:n
    l₁ = Line(vother[i], vother[mod1(i+1, n)])

    m = length(v)
    newv = Point{Dim,T}[]
    
    for j in 1:m
      p₁, p₂ = v[j], v[mod1(j+1, m)]
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

  isempty(v) ? nothing : PolyArea(Ring(v))
end
