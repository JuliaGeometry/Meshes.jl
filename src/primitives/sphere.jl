# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Sphere(center, radius)

A sphere with `center` and `radius`.

See also [`Ball`](@ref).
"""
struct Sphere{Dim,T} <: Primitive{Dim,T}
  center::Point{Dim,T}
  radius::T
end

Sphere(center::Point{Dim,T}, radius) where {Dim,T} = Sphere(center, T(radius))

Sphere(center::Tuple, radius) = Sphere(Point(center), radius)

"""
    Sphere(p1, p2, p3)

A 2D sphere passing through points `p1`, `p2` and `p3`.
"""
function Sphere(p1::Point{2}, p2::Point{2}, p3::Point{2})
  x1, y1 = p2 - p1
  x2, y2 = p3 - p2
  c1 = centroid(Segment(p1, p2))
  c2 = centroid(Segment(p2, p3))
  l1 = Line(c1, c1 + Vec(y1, -x1))
  l2 = Line(c2, c2 + Vec(y2, -x2))
  center = l1 ∩ l2
  radius = norm(center - p2)
  Sphere(center, radius)
end

Sphere(p1::Tuple, p2::Tuple, p3::Tuple) = Sphere(Point(p1), Point(p2), Point(p3))

"""
    Sphere(p1, p2, p3, p4)

A 3D sphere passing through points `p1`, `p2`, `p3` and `p4`.
"""
function Sphere(p1::Point{3}, p2::Point{3}, p3::Point{3}, p4::Point{3})
  x1 = coordinates(p1)
  x2 = coordinates(p2)
  x3 = coordinates(p3)
  x4 = coordinates(p4)
  X = hcat(x1, x2, x3, x4)
  a = 2 * det(vcat(X, ones((1, 4))))
  T = typeof(a)
  if isapprox(a, zero(T); atol = atol(T))
    error("The four points are coplanar.")
  end
  q = [x1 ⋅ x1, x2 ⋅ x2, x3 ⋅ x3, x4 ⋅ x4] / a
  y1 = X[1, :]
  y2 = X[2, :]
  y3 = X[3, :]
  d1 = det(hcat(q, y2, y3, ones(4)))
  d2 = -det(hcat(q, y1, y3, ones(4)))
  d3 = det(hcat(q, y1, y2, ones(4)))
  c = det(hcat(q, y1, y2, y3))
  radius = sqrt(d1^2 + d2^2 + d3^2 - 2 * c) 
  center = Point(d1, d2, d3)
  Sphere(center, radius)
end

Sphere(p1::Tuple, p2::Tuple, p3::Tuple, p4::Tuple) = 
  Sphere(Point(p1), Point(p2), Point(p3), Point(p4))

paramdim(::Type{<:Sphere{Dim}}) where {Dim} = Dim - 1

isconvex(::Type{<:Sphere}) = false

isperiodic(::Type{<:Sphere{Dim}}) where {Dim} = ntuple(i->true, Dim - 1)

center(s::Sphere) = s.center

radius(s::Sphere) = s.radius

# https://en.wikipedia.org/wiki/N-sphere#Volume_and_surface_area
function measure(s::Sphere{Dim}) where {Dim}
  r, n = s.radius, Dim
  2π^(n/2)*r^(n-1) / gamma(n/2)
end

Base.length(s::Sphere{2}) = measure(s)

area(s::Sphere{3}) = measure(s)

boundary(::Sphere) = nothing

perimeter(::Sphere{Dim,T}) where {Dim,T} = zero(T)

function Base.in(p::Point, s::Sphere)
  x = coordinates(p)
  c = coordinates(s.center)
  r = s.radius
  sum(abs2, x - c) ≈ r^2
end
