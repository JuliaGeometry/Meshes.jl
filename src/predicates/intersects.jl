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

intersects(r::Ray, b::Box) = !isnothing(r ∩ b)

intersects(b::Box, r::Ray) = intersects(r, b)

intersects(r::Ray, t::Triangle) = !isnothing(r ∩ t)

intersects(t::Triangle, r::Ray) = intersects(r, t)

function intersects(r::Ray, s::Sphere)
  u = center(s) - r(0)
  h = norm(u)
  radius(s) > h && return true
  v = r(1) - r(0)
  abs(∠(u, v)) < asin(radius(s) / h)
end

intersects(s::Sphere, r::Ray) = intersects(r, s)

intersects(r::Ray, b::Ball) = intersects(r, boundary(b))

intersects(b::Ball, r::Ray) = intersects(r, b)

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
    else
      d = gjk!(Val(Dim), O, points)
      isnothing(d) && return true
    end
  end
end

"""
    gjk!(::Val{Dim}, origin, points) where Dim

Perform one iteration of the GJK algorithm.

It returns `nothing` if the `Dim`-simplex represented by `points`
contains the origin point. Otherwise, it returns a vector with
the direction for searching the next point.

If the simplex is complete, it removes one point from the set to
make room for the next point. A complete simplex must have `Dim + 1` points. 

See also [`intersects`](@ref).
"""
function gjk! end

function gjk!(::Val{2}, origin, points)
  @assert length(points) == 3
  # triangle simplex case
  C, B, A = points
  AB = B - A
  AC = C - A
  AO = origin - A
  ABᵀ = AC × AB × AB
  ACᵀ = AB × AC × AC
  T = coordtype(origin)
  if ABᵀ ⋅ AO > zero(T)
    popat!(points, 1) # pop C
    return ABᵀ
  elseif ACᵀ ⋅ AO > zero(T)
    popat!(points, 2) # pop B
    return ACᵀ
  else
    return nothing
  end
end

function gjk!(::Val{3}, origin, points)
  if length(points) == 3
    # triangle case
    C, B, A = points
    AB = B - A
    AC = C - A
    AO = origin - A
    ABCᵀ = AB × AC
    if ABCᵀ ⋅ AO < 0
      points[1], points[2] = points[2], points[1]
      ABCᵀ = -ABCᵀ
    end
    return ABCᵀ
  else
    # tetrahedron simplex case
    #      A
    #    / | \
    #   /  D  \
    #  / /   \ \
    # B ------- C
    # Simplex faces (with normal vectors pointing to the centroid):
    # ABC, ADB, BDC, ACD
    # (AXY = AX × AY)
    # ABC normal vector points to vertex D
    D, C, B, A = points
    AB = B - A
    AC = C - A
    AD = D - A
    AO = origin - A
    ABCᵀ = AB × AC
    ADBᵀ = AD × AB
    ACDᵀ = AC × AD
    T = coordtype(origin)
    if ABCᵀ ⋅ AO > zero(T)
      popat!(points, 1) # pop D
      return ABCᵀ
    elseif ADBᵀ ⋅ AO > zero(T)
      popat!(points, 2) # pop C
      return ADBᵀ
    elseif ACDᵀ ⋅ AO > zero(T)
      popat!(points, 3) # pop B
      return ACDᵀ
    else
      return nothing
    end
  end
end

intersects(m::Multi, g::Geometry) = intersects(parent(m), [g])

intersects(g::Geometry, m::Multi) = intersects(m, g)

intersects(m₁::Multi, m₂::Multi) = intersects(parent(m₁), parent(m₂))

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

intersects(c::Chain, m::Multi) = intersects(segments(c), parent(m))

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
