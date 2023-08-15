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

  for s1 in segments(other)
    n = length(v)
    newv = []

    for i in 1:n
      v₁, v₂ = v[i], v[i%n+1]
      s2 = Segment(v₁, v₂)

      # assuming convex clockwise other
      isv₁inside = (sideof(v₁, s1) != :LEFT)
      isv₂inside = (sideof(v₂, s1) != :LEFT)

      if isv₁inside && isv₂inside
        push!(newv, v₁)
      elseif isv₁inside && !isv₂inside
        push!(newv, v₁)
        push!(newv, s1 ∩ s2)
      elseif !isv₁inside && isv₂inside
        push!(newv, s1 ∩ s2)
      end
    end
    v = newv
  end
  v
end
