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

function clip(poly, other, algo) 
  r = [clip(ring, boundary(other), algo) for ring in rings(poly)]
  # TODO: filter nothing rings
  if hasholes(poly)
    PolyArea(r[1], r[2:end])
  else
    PolyArea(r[1])
  end
end

function clip(poly::Ring{Dim,T}, other::Ring{Dim,T}, ::SutherlandHodgman) where {Dim,T}
  o = vertices(other)
  n = length(o)
  v = vertices(poly)

  for i in 1:n
    l₁ = Line(o[i], o[mod1(i+1, n)])

    m = length(v)
    u = Point{Dim,T}[]
    
    for j in 1:m
      p₁, p₂ = v[j], v[mod1(j+1, m)]
      l₂ = Line(p₁, p₂)

      # assuming convex clockwise other
      isinside₁ = (sideof(p₁, l₁) != :LEFT)
      isinside₂ = (sideof(p₂, l₁) != :LEFT)

      if isinside₁ && isinside₂
        push!(u, p₁)
      elseif isinside₁ && !isinside₂
        push!(u, p₁)
        push!(u, l₁ ∩ l₂)
      elseif !isinside₁ && isinside₂
        push!(u, l₁ ∩ l₂)
      end
    end
    v = u
  end

  isempty(v) ? nothing : Ring(v)
end
