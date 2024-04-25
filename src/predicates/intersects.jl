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

function intersects(g₁::Geometry{Dim}, g₂::Geometry{Dim}) where {Dim}
  𝒬 = coordtype(g₁)

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
  d = c₁ ≈ c₂ ? rand(Vec{Dim,𝒬}) : c₂ - c₁

  # first point in Minkowski difference
  P = minkowskipoint(g₁, g₂, d)

  # origin of coordinate system
  O = minkowskiorigin(Dim, 𝒬)

  # initialize simplex vertices
  points = [P]

  # move towards the origin
  d = O - P
  while true
    P = minkowskipoint(g₁, g₂, d)
    if (P - O) ⋅ d < zero(𝒬)
      return false
    end
    push!(points, P)

    d = gjk!(O, points)
    isnothing(d) && return true
  end
end

"""
    gjk!(O::Point{Dim,T}, points) where {Dim,T}

Perform one iteration of the GJK algorithm.

It returns `nothing` if the `Dim`-simplex represented by `points`
contains the origin point `O`. Otherwise, it returns a vector with
the direction for searching the next point.

If the simplex is complete, it removes one point from the set to
make room for the next point. A complete simplex must have `Dim + 1` points. 

See also [`intersects`](@ref).
"""
function gjk! end

function gjk!(O::Point{2}, points)
  𝒬 = coordtype(O)
  # line segment case
  if length(points) == 2
    B, A = points
    AB = B - A
    AO = O - A
    d = perpendicular(AB, AO)
  else
    # triangle simplex case
    C, B, A = points
    AB = B - A
    AC = C - A
    AO = O - A
    ABᵀ = -perpendicular(AB, AC)
    ACᵀ = -perpendicular(AC, AB)
    if ABᵀ ⋅ AO > zero(𝒬)
      popat!(points, 1) # pop C
      d = ABᵀ
    elseif ACᵀ ⋅ AO > zero(𝒬)
      popat!(points, 2) # pop B
      d = ACᵀ
    else
      d = nothing
    end
  end
  d
end

function gjk!(O::Point{3}, points)
  𝒬 = coordtype(O)
  # line segment case
  if length(points) == 2
    B, A = points
    AB = B - A
    AO = O - A
    d = perpendicular(AB, AO)
  elseif length(points) == 3
    # triangle case
    C, B, A = points
    AB = B - A
    AC = C - A
    AO = O - A
    ABCᵀ = AB × AC
    if ABCᵀ ⋅ AO < 0
      points[1], points[2] = points[2], points[1]
      ABCᵀ = -ABCᵀ
    end
    d = ABCᵀ
  else
    # tetrahedron simplex case
    #      A
    #    / | \
    #   /  D  \
    #  / /   \ \
    # C ------- B
    # Simplex faces (with normal vectors pointing away from the centroid):
    # ABC, ADB, BDC, ACD
    # (AXY = AX × AY)
    # ACB normal vector points to vertex D
    # ABC normal vector points in the opposite direction
    D, C, B, A = points
    AB = B - A
    AC = C - A
    AD = D - A
    AO = O - A
    ABCᵀ = AB × AC
    ADBᵀ = AD × AB
    ACDᵀ = AC × AD
    if ABCᵀ ⋅ AO > zero(𝒬)
      popat!(points, 1) # pop D
      d = ABCᵀ
    elseif ADBᵀ ⋅ AO > zero(𝒬)
      popat!(points, 2) # pop C
      d = ADBᵀ
    elseif ACDᵀ ⋅ AO > zero(𝒬)
      popat!(points, 3) # pop B
      d = ACDᵀ
    else
      d = nothing
    end
  end
  d
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
minkowskipoint(g₁::Geometry, g₂::Geometry, d) = Point(supportfun(g₁, d) - supportfun(g₂, -d))

# origin of coordinate system
minkowskiorigin(Dim, T) = Point(ntuple(i -> zero(T), Dim))

# find a vector perpendicular to `v` using vector `d` as some direction hint
# expect that `perpendicular(v, d) ⋅ d ≥ 0` or, in other words,
# that the angle between the result vector and `d` is less or equal than 90º
function perpendicular(v::Vec{2,L}, d::Vec{2,L}) where {L}
  a = Vec(v[1], v[2], zero(L))
  b = Vec(d[1], d[2], zero(L))
  r = a × b × a
  Vec(r[1], r[2])
end

perpendicular(v::Vec{3}, d::Vec{3}) = v × d × v
