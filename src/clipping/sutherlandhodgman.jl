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
  # make sure other ring is CCW
  occw = orientation(other) == CCW ? other : reverse(other)

  v = vertices(ring)
  o = vertices(occw)

  for i in 1:length(o)
    lₒ = Line(o[i], o[i + 1])

    n = length(v)

    u = Point{Dim,T}[]
    for j in 1:n
      v₁ = v[j]
      v₂ = v[mod1(j + 1, n)]
      lᵣ = Line(v₁, v₂)

      isinside₁ = (sideof(v₁, lₒ) != RIGHT)
      isinside₂ = (sideof(v₂, lₒ) != RIGHT)

      if isinside₁ && isinside₂
        push!(u, v₁)
      elseif isinside₁ && !isinside₂
        push!(u, v₁)
        push!(u, lₒ ∩ lᵣ)
      elseif !isinside₁ && isinside₂
        push!(u, lₒ ∩ lᵣ)
      end
    end

    v = u
  end

  isempty(v) ? nothing : Ring(unique(v))
end
