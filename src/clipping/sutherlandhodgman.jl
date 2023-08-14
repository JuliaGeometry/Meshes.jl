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


clip(poly::T, window::Polygon, ::SutherlandHodgman) where {T <: Polygon} = T(clip(vertices(poly), segments(window)))

function clip(poly::PolyArea, window::Polygon, ::SutherlandHodgman)
  r = map(rings(poly)) do ring
    clip(ring, window, ::SutherlandHodgman)
  end
  PolyArea(first(r), r[2:end])
end

# ---------------
# IMPLEMENTATION
# ---------------

function clip(v::AbstractVector{Point}, window::AbstractVector{Segment}, ::SutherlandHodgman)
  # clip each segment of the window one per time
  for s in window
    v = clip(v, s, ::SutherlandHodgman)
  end
  v
end

function clip(v::AbstractVector{Point}, window::Segment, ::SutherlandHodgman)
  n = length(v)
  new_v = []

  for i in 1:n
    p1 = v[i]
    p2 = v[(i%n)+1]

    # assuming convex clockwise window
    p1_visible = (sideof(p1, window) != :LEFT) 
    p2_visible = (sideof(p2, window) != :LEFT)

    if p1_visible && p2_visible
      push!(new_v, p1)
    elseif p1_visible && !p2_visible
      p_intersection = Segment(p1, p2) ∩ window
      push!(new_v, p1)
      push!(new_v, p_intersection)
    elseif !p1_visible && p2_visible
      p_intersection = Segment(p1, p2) ∩ window
      push!(new_v, p_intersection)
    end
  end
  new_v
end
