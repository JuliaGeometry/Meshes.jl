# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Sphere(center, radius)

A sphere with `center` and `radius`.

See also [`Ball`](@ref).
"""
struct Sphere{M<:Manifold,C<:CRS,ℒ<:Len} <: Primitive{M,C}
  center::Point{M,C}
  radius::ℒ
  Sphere(center::Point{M,C}, radius::ℒ) where {M<:Manifold,C<:CRS,ℒ<:Len} = new{M,C,float(ℒ)}(center, radius)
end

Sphere(center::Point, radius) = Sphere(center, addunit(radius, u"m"))

Sphere(center::Tuple, radius) = Sphere(Point(center), radius)

Sphere(center::Point) = Sphere(center, oneunit(lentype(center)))

Sphere(center::Tuple) = Sphere(Point(center))

"""
    Sphere(p1, p2, p3)

A 2D sphere passing through points `p1`, `p2` and `p3`.
"""
function Sphere(p1::Point, p2::Point, p3::Point)
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
function Sphere(p1::Point, p2::Point, p3::Point, p4::Point)
  v1 = p1 - p4
  v2 = p2 - p4
  v3 = p3 - p4
  V = volume(Tetrahedron(p1, p2, p3, p4))
  r⃗ = ((v3 ⋅ v3) * (v1 × v2) + (v2 ⋅ v2) * (v3 × v1) + (v1 ⋅ v1) * (v2 × v3)) / 12V
  center = p4 + Vec(r⃗)
  radius = norm(r⃗)
  Sphere(center, radius)
end

Sphere(p1::Tuple, p2::Tuple, p3::Tuple, p4::Tuple) = Sphere(Point(p1), Point(p2), Point(p3), Point(p4))

paramdim(::Type{<:Sphere{𝔼{Dim}}}) where {Dim} = Dim - 1

paramdim(::Type{<:Sphere{🌐}}) = 1

center(s::Sphere) = s.center

radius(s::Sphere) = s.radius

==(s₁::Sphere, s₂::Sphere) = s₁.center == s₂.center && s₁.radius == s₂.radius

Base.isapprox(s₁::Sphere, s₂::Sphere; atol=atol(lentype(s₁)), kwargs...) =
  isapprox(s₁.center, s₂.center; atol, kwargs...) && isapprox(s₁.radius, s₂.radius; atol, kwargs...)

(s::Sphere)(uv...) = Ball(center(s), radius(s))(1, uv...)
