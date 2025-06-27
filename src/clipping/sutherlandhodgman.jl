# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SutherlandHodgmanClipping()

The Sutherland-Hodgman algorithm for clipping polygons.

## References

* Sutherland, I.E. & Hodgman, G.W. 1974. [Reentrant Polygon Clipping]
  (https://dl.acm.org/doi/pdf/10.1145/360767.360802)

### Notes

The algorithm assumes that the clipping geometry is convex.
"""
struct SutherlandHodgmanClipping <: ClippingMethod end

function clip(poly::Polygon, other::Geometry, method::SutherlandHodgmanClipping)
  c = [clip(ring, boundary(other), method) for ring in rings(poly)]
  r = [r for r in c if !isnothing(r)]
  isempty(r) ? nothing : PolyArea(r)
end

function clip(ring::Ring, other::Ring, ::SutherlandHodgmanClipping)
  # make sure other ring is CCW
  occw = orientation(other) == CCW ? other : reverse(other)

  # vertices as circular vectors
  vᵣ = vertices(ring)
  vₒ = vertices(occw)

  for j in eachindex(vₒ)
    lₒ = Line(vₒ[j], vₒ[j + 1])

    # retain vertices from vᵣ that satisfy
    # Sutherland-Hodgeman criteria
    p = empty(vᵣ)
    for i in eachindex(vᵣ)
      p₁ = vᵣ[i]
      p₂ = vᵣ[i + 1]
      lᵣ = Line(p₁, p₂)

      isinside₁ = (sideof(p₁, lₒ) != RIGHT)
      isinside₂ = (sideof(p₂, lₒ) != RIGHT)

      if isinside₁ && isinside₂
        push!(p, p₁)
      elseif isinside₁ && !isinside₂
        push!(p, p₁)
        push!(p, intersectpoint(lᵣ, lₒ))
      elseif !isinside₁ && isinside₂
        push!(p, intersectpoint(lᵣ, lₒ))
      end
    end

    # update list of vertices and continue
    vᵣ = p
  end

  # return appropriate object
  isempty(vᵣ) ? nothing : Ring(unique(vᵣ))
end

# helper function to find any intersection point
# between crossing or overlapping lines
function intersectpoint(l₁::Line, l₂::Line)
  λ(I) = type(I) == Overlapping ? l₁(0) : get(I)
  intersection(λ, l₁, l₂)
end
