# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SutherlandHodgmanClipping()

The Sutherland-Hodgman algorithm for clipping polygons.

## References

* Sutherland, I.E. & Hodgman, G.W. 1974. [Reentrant Polygon
  Clipping](https://dl.acm.org/doi/pdf/10.1145/360767.360802)

### Notes

The algorithm assumes that the other geometry [`isconvex`](@ref).
"""
struct SutherlandHodgmanClipping <: ClippingMethod end

function clip(subject::Polygon, other::Polygon, ::SutherlandHodgmanClipping)
  orings = rings(other)
  if length(orings) > 1
    throw(ArgumentError("Sutherland-Hodgman requires convex clipping polygon"))
  end
  srings = rings(subject)
  crings = empty(srings)
  for sring in srings
    cring = _shclip(sring, first(orings))
    isnothing(cring) || push!(crings, cring)
  end
  isempty(crings) ? nothing : PolyArea(crings)
end

function _shclip(ring::Ring, other::Ring)
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
        push!(p, _shpoint(lᵣ, lₒ))
      elseif !isinside₁ && isinside₂
        push!(p, _shpoint(lᵣ, lₒ))
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
function _shpoint(l₁::Line, l₂::Line)
  λ(I) = type(I) == Overlapping ? l₁(0) : get(I)
  intersection(λ, l₁, l₂)
end
