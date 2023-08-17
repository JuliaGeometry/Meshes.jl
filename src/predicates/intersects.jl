# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    intersects(geometry₁, geometry₂)

Tells whether or not `geometry₁` and `geometry₂` intersect.

## References

* Gilbert, E., Johnson, D., Keerthi, S. 1988. [A fast
  Procedure for Computing the Distance Between Complex
  Objects in Three-Dimensional Space]
  (https://ieeexplore.ieee.org/document/2083)

### Notes

* The fallback algorithm works with any geometry that has
  a well-defined [`supportfun`](@ref).
"""
function intersects end

intersects(g) = Base.Fix2(intersects, g)

intersects(p₁::Point, p₂::Point) = p₁ == p₂

intersects(s₁::Segment, s₂::Segment) = !isnothing(s₁ ∩ s₂)

intersects(b₁::Box, b₂::Box) = !isnothing(b₁ ∩ b₂)

intersects(r::Ray, t::Triangle) = !isnothing(r ∩ t)

intersects(t::Triangle, r::Ray) = intersects(r, t)

intersects(p::Point, g::Geometry) = p ∈ g

intersects(g::Geometry, p::Point) = intersects(p, g)

intersects(c::Chain, s::Segment) = intersects(segments(c), [s])

intersects(s::Segment, c::Chain) = intersects(c, s)

intersects(c₁::Chain, c₂::Chain) = intersects(segments(c₁), segments(c₂))

intersects(c::Chain, g::Geometry) = any(∈(g), vertices(c)) || intersects(c, boundary(g))

intersects(g::Geometry, c::Chain) = intersects(c, g)

function intersects(g₁::Geometry{Dim,T}, g₂::Geometry{Dim,T}) where {Dim,T}
  # must have intersection of bounding boxes
  intersects(boundingbox(g₁), boundingbox(g₂)) || return false

  # handle non-convex geometries
  if !isconvex(g₁)
    d₁ = simplexify(g₁)
    return intersects(d₁, g₂)
  elseif !isconvex(g₂)
    d₂ = simplexify(g₂)
    return intersects(g₁, d₂)
  end

  # initial direction
  c₁, c₂ = centroid(g₁), centroid(g₂)
  d = c₁ ≈ c₂ ? rand(Vec{Dim,T}) : c₂ - c₁

  # first point in Minkowski difference
  P = minkowskipoint(g₁, g₂, d)

  # origin of coordinate system
  O = minkowskiorigin(Dim, T)

  # initialize simplex vertices
  points = [P]

  # move towards the origin
  d = O - P
  while true
    P = minkowskipoint(g₁, g₂, d)
    if (P - O) ⋅ d < zero(T)
      return false
    end
    push!(points, P)

    # line segment case
    if length(points) == 2
      B, A = points
      AB = B - A
      AO = O - A
      d = AB × AO × AB
    else # simplex case
      C, B, A = points
      AB = B - A
      AC = C - A
      AO = O - A
      ABᵀ = AC × AB × AB
      ACᵀ = AB × AC × AC
      if ABᵀ ⋅ AO > zero(T)
        popat!(points, 1) # pop C
        d = ABᵀ
      elseif ACᵀ ⋅ AO > zero(T)
        popat!(points, 2) # pop B
        d = ACᵀ
      else
        return true
      end
    end
  end
end

intersects(m::Multi, g::Geometry) = intersects(collect(m), [g])

intersects(g::Geometry, m::Multi) = intersects(m, g)

intersects(m₁::Multi, m₂::Multi) = intersects(collect(m₁), collect(m₂))

intersects(d::Domain, g::Geometry) = intersects(d, [g])

intersects(g::Geometry, d::Domain) = intersects(d, g)

# fallback with iterators of geometries
function intersects(geoms₁, geoms₂)
  for g₁ in geoms₁, g₂ in geoms₂
    intersects(g₁, g₂) && return true
  end
  return false
end

# -------------------------
# solve method ambiguities
# -------------------------

intersects(p::Point, c::Chain) = p ∈ c

intersects(c::Chain, p::Point) = intersects(p, c)

intersects(p::Point, m::Multi) = p ∈ m

intersects(m::Multi, p::Point) = intersects(p, m)

intersects(c::Chain, m::Multi) = intersects(segments(c), collect(m))

intersects(m::Multi, c::Chain) = intersects(c, m)

# ------------------
# utility functions
# ------------------

# support point in Minkowski difference
function minkowskipoint(g₁::Geometry{Dim,T}, g₂::Geometry{Dim,T}, d) where {Dim,T}
  n = Vec{Dim,T}(d[1:Dim])
  v = supportfun(g₁, n) - supportfun(g₂, -n)
  Point(ntuple(i -> i ≤ Dim ? v[i] : zero(T), max(Dim, 3)))
end

# origin of coordinate system
minkowskiorigin(Dim, T) = Point(ntuple(i -> zero(T), max(Dim, 3)))
