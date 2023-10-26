# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SutherlandHodgman()

The Sutherland-Hodgman algorithm for clipping polygons.

## References

* Sutherland, I.E. & Hodgman, G.W. 1974. [Reentrant Polygon Clipping]
  (https://dl.acm.org/doi/pdf/10.1145/360767.360802)

### Notes

* The algorithm assumes that the clipping geometry is convex.
"""
struct SutherlandHodgman <: ClippingMethod end

function clip(poly::Polygon, other::Geometry, method::SutherlandHodgman)
  c = [clip(ring, boundary(other), method) for ring in rings(poly)]
  r = [r for r in c if !isnothing(r)]
  isempty(r) ? nothing : PolyArea(r)
end

function clip(ring::Ring{Dim,T}, other::Ring{Dim,T}, ::SutherlandHodgman) where {Dim,T}
  v = vertices(ring)

  o = orientation(other) == CCW ? vertices(other) : reverse(vertices(other))
  n = length(o)

  for i in 1:n
    lₒ = Line(o[i], o[mod1(i + 1, n)])

    m = length(v)
    u = Point{Dim,T}[]

    for j in 1:m
      p₁ = v[j]
      p₂ = v[mod1(j + 1, m)]
      lᵣ = Line(p₁, p₂)

      isinside₁ = (sideof(p₁, lₒ) != RIGHT)
      isinside₂ = (sideof(p₂, lₒ) != RIGHT)

      if isinside₁ && isinside₂
        push!(u, p₁)
      elseif isinside₁ && !isinside₂
        push!(u, p₁)
        push!(u, lₒ ∩ lᵣ)
      elseif !isinside₁ && isinside₂
        push!(u, lₒ ∩ lᵣ)
      end
    end

    v = u
  end

  isempty(v) ? nothing : Ring(unique(v))
end
