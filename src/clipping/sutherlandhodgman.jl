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

function clip(ring::Ring, other::Ring, ::SutherlandHodgman)
  # make sure other ring is CCW
  occw = orientation(other) == CCW ? other : reverse(other)

  r = vertices(ring)
  o = vertices(occw)

  for i in 1:length(o)
    lₒ = Line(o[i], o[i + 1])

    n = length(r)

    u = Vector{eltype(r)}()
    for j in 1:n
      r₁ = r[j]
      r₂ = r[mod1(j + 1, n)]
      lᵣ = Line(r₁, r₂)

      isinside₁ = (sideof(r₁, lₒ) != RIGHT)
      isinside₂ = (sideof(r₂, lₒ) != RIGHT)

      if isinside₁ && isinside₂
        push!(u, r₁)
      elseif isinside₁ && !isinside₂
        push!(u, r₁)
        push!(u, intersectpoint(lᵣ, lₒ))
      elseif !isinside₁ && isinside₂
        push!(u, intersectpoint(lᵣ, lₒ))
      end
    end

    r = u
  end

  isempty(r) ? nothing : Ring(unique(r))
end

# helper function to find any intersection point
# between crossing or overlapping lines
function intersectpoint(l₁::Line, l₂::Line)
  intersection(l₁, l₂) do I
    if type(I) == Crossing # get intersection point
      get(I)
    elseif type(I) == Overlapping # perturb line and retry
      ℒ = lentype(l₁)
      T = numtype(ℒ)
      δx = rand((1, -1)) * rand(T) * atol(ℒ)
      δy = rand((1, -1)) * rand(T) * atol(ℒ)
      l₁′ = Line(l₁(0), l₁(1) + Vec(δx, δy))
      intersectpoint(l₁′, l₂)
    end
  end
end
