# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Sphere(center, radius)

A sphere with `center` and `radius`.

See also [`Ball`](@ref).
"""
struct Sphere{Dim,P<:Point{Dim},L<:Len} <: Primitive{Dim}
  center::P
  radius::L
  Sphere(center::P, radius::L) where {Dim,P<:Point{Dim},L<:Len} = new{Dim,P,float(L)}(center, radius)
end

Sphere(center::Point, radius) = Sphere(center, addunit(radius, u"m"))

Sphere(center::Tuple, radius) = Sphere(Point(center), radius)

Sphere(center::Point) = Sphere(center, 1.0)

Sphere(center::Tuple) = Sphere(Point(center))

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
  v1 = p1 - p4
  v2 = p2 - p4
  v3 = p3 - p4
  V = volume(Tetrahedron(p1, p2, p3, p4))
  r⃗ = ((v3 ⋅ v3) * (v1 × v2) + (v2 ⋅ v2) * (v3 × v1) + (v1 ⋅ v1) * (v2 × v3)) / 12V
  center = p4 + r⃗
  radius = norm(r⃗)
  Sphere(center, radius)
end

Sphere(p1::Tuple, p2::Tuple, p3::Tuple, p4::Tuple) = Sphere(Point(p1), Point(p2), Point(p3), Point(p4))

paramdim(::Type{<:Sphere{Dim}}) where {Dim} = Dim - 1

coordtype(::Type{<:Sphere{Dim,P}}) where {Dim,P} = coordtype(P)

center(s::Sphere) = s.center

radius(s::Sphere) = s.radius

function (s::Sphere{2,P,L})(φ) where {P,L}
  if (φ < 0 || φ > 1)
    throw(DomainError(φ, "s(φ) is not defined for φ outside [0, 1]."))
  end
  c = s.center
  r = s.radius
  sφ, cφ = sincospi(2 * L(φ))
  x = r * cφ
  y = r * sφ
  c + Vec(x, y)
end

function (s::Sphere{3,P,L})(θ, φ) where {P,L}
  if (θ < 0 || θ > 1) || (φ < 0 || φ > 1)
    throw(DomainError((θ, φ), "s(θ, φ) is not defined for θ, φ outside [0, 1]²."))
  end
  c = s.center
  r = s.radius
  sθ, cθ = sincospi(L(θ))
  sφ, cφ = sincospi(2 * L(φ))
  x = r * sθ * cφ
  y = r * sθ * sφ
  z = r * cθ
  c + Vec(x, y, z)
end

# TODO
# Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Sphere{Dim,T}}) where {Dim,T} =
#   Sphere(rand(rng, Point{Dim,T}), rand(rng, T))
