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

The fallback algorithm works with any geometry that has
a well-defined [`supportfun`](@ref).
"""
function intersects end

intersects(i::Intersection) = type(i) !== NotIntersecting

intersects(g) = Base.Fix2(intersects, g)

intersects(p₁::Point, p₂::Point) = p₁ == p₂

intersects(s₁::Segment, s₂::Segment) = intersects(intersection(s₁, s₂))

intersects(b₁::Box, b₂::Box) = intersects(intersection(b₁, b₂))

intersects(r::Ray, b::Box) = intersects(intersection(r, b))

intersects(b::Box, r::Ray) = intersects(r, b)

intersects(r::Ray, t::Triangle) = intersects(intersection(r, t))

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

intersects(c::Chain, g::Geometry) = any(∈(g), eachvertex(c)) || intersects(c, boundary(g))

intersects(g::Geometry, c::Chain) = intersects(c, g)

function intersects(g₁::Geometry, g₂::Geometry)
  Dim = embeddim(g₁)
  ℒ = lentype(g₁)

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
  d = c₁ ≈ c₂ ? rand(Vec{Dim,ℒ}) : c₂ - c₁

  # first point in Minkowski difference
  P = minkowskipoint(g₁, g₂, d)

  # origin of coordinate system
  O = minkowskiorigin(Dim, ℒ)

  # initialize simplex vertices
  points = [P]

  # move towards the origin
  d = O - P
  while true
    P = minkowskipoint(g₁, g₂, d)
    if isnegative((P - O) ⋅ d)
      return false
    end
    push!(points, P)

    d = gjk!(O, points)
    isnothing(d) && return true
  end
end

"""
    gjk!(O::Point{Dim}, points) where {Dim}

Perform one iteration of the GJK algorithm.

It returns `nothing` if the `Dim`-simplex represented by `points`
contains the origin point `O`. Otherwise, it returns a vector with
the direction for searching the next point.

If the simplex is complete, it removes one point from the set to
make room for the next point. A complete simplex must have `Dim + 1` points.

See also [`intersects`](@ref).
"""
gjk!(O::Point, points) = _gjk!(Val(embeddim(O)), O, points)

function _gjk!(::Val{2}, O, points)
  # line segment case
  if length(points) == 2
    B, A = points
    AB = B - A
    AO = O - A
    d = perphint(AB, AO)
  else
    # triangle simplex case
    C, B, A = points
    AB = B - A
    AC = C - A
    AO = O - A
    ABᵀ = -perphint(AB, AC)
    ACᵀ = -perphint(AC, AB)
    if ispositive(ABᵀ ⋅ AO)
      popat!(points, 1) # pop C
      d = ABᵀ
    elseif ispositive(ACᵀ ⋅ AO)
      popat!(points, 2) # pop B
      d = ACᵀ
    else
      d = nothing
    end
  end
  d
end

function _gjk!(::Val{3}, O, points)
  # line segment case
  if length(points) == 2
    B, A = points
    AB = B - A
    AO = O - A
    d = perphint(AB, AO)
  elseif length(points) == 3
    # triangle case
    C, B, A = points
    AB = B - A
    AC = C - A
    AO = O - A
    ABCᵀ = ucross(AB, AC)
    if isnegative(ABCᵀ ⋅ AO)
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
    ABCᵀ = ucross(AB, AC)
    ADBᵀ = ucross(AD, AB)
    ACDᵀ = ucross(AC, AD)
    if ispositive(ABCᵀ ⋅ AO)
      popat!(points, 1) # pop D
      d = ABCᵀ
    elseif ispositive(ADBᵀ ⋅ AO)
      popat!(points, 2) # pop C
      points[1], points[2] = points[2], points[1]
      d = ADBᵀ
    elseif ispositive(ACDᵀ ⋅ AO)
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
minkowskipoint(g₁::Geometry, g₂::Geometry, d) = withcrs(g₁, supportfun(g₁, d) - supportfun(g₂, -d))

# origin of coordinate system
minkowskiorigin(Dim, ℒ) = Point(ntuple(i -> zero(ℒ), Dim))

# find a vector perpendicular to `v` using vector `d` as some direction hint
# expect that `perphint(v, d) ⋅ d ≥ 0` or, in other words,
# that the angle between the result vector and `d` is less or equal than 90º
function perphint(v::Vec{2,ℒ}, d::Vec{2,ℒ}) where {ℒ}
  a = Vec(v[1], v[2], zero(ℒ))
  b = Vec(d[1], d[2], zero(ℒ))
  r = ucross(a, b, a)
  Vec(r[1], r[2])
end

perphint(v::Vec{3,ℒ}, d::Vec{3,ℒ}) where {ℒ} = ucross(v, d, v)

