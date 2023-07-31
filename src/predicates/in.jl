# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    point ∈ geometry

Tells whether or not the `point` is in the `geometry`.
"""
function Base.in(::Point, ::Geometry) end

Base.in(p₁::Point, p₂::Point) = p₁ == p₂

function Base.in(p::Point{Dim,T}, s::Segment{Dim,T}) where {Dim,T}
  # given collinear points (a, b, p), the point p intersects
  # segment ab if and only if vectors satisfy 0 ≤ ap ⋅ ab ≤ ||ab||²
  a, b = s.vertices
  ab, ap = b - a, p - a
  iscollinear(a, b, p) && zero(T) ≤ ab ⋅ ap ≤ ab ⋅ ab
end

Base.in(p::Point, r::Ray) = p ∈ Line(r.p, r.p + r.v) && (p - r.p) ⋅ r.v ≥ 0

function Base.in(p::Point, l::Line)
  w = norm(l.b - l.a)
  d = evaluate(Euclidean(), p, l)
  # d ≈ 0.0 will be too precise, and d < atol{T} can't scale.
  d + w ≈ w
end

Base.in(p::Point, c::Chain) = any(s -> p ∈ s, segments(c))

Base.in(pt::Point{3,T}, pl::Plane{T}) where {T} = isapprox(normal(pl) ⋅ (pt - pl(0, 0)), zero(T), atol=atol(T))

function Base.in(p::Point{Dim}, b::Box{Dim}) where {Dim}
  l, u = coordinates.((b.min, b.max))
  x = coordinates(p)
  for i in 1:Dim
    l[i] ≤ x[i] && x[i] ≤ u[i] || return false
  end
  true
end

function Base.in(p::Point{Dim,T}, b::Ball{Dim,T}) where {Dim,T}
  c = b.center
  r = b.radius
  s = norm(p - c)
  s < r || isapprox(s, r, atol=atol(T))
end

function Base.in(p::Point, s::Sphere)
  x = coordinates(p)
  c = coordinates(s.center)
  r = s.radius
  sum(abs2, x - c) ≈ r^2
end

function Base.in(p::Point, d::Disk)
  p ∉ d.plane && return false
  s² = sum(abs2, p - center(d))
  r² = radius(d)^2
  s² ≤ r²
end

function Base.in(p::Point{3,T}, c::Circle{T}) where {T}
  p ∉ c.plane && return false
  s² = sum(abs2, p - center(c))
  r² = radius(c)^2
  isapprox(s², r², atol=atol(T)^2)
end

function Base.in(p::Point{3}, c::Cylinder)
  b = c.bot(0, 0)
  t = c.top(0, 0)
  a = t - b
  (p - b) ⋅ a ≥ 0 || return false
  (p - t) ⋅ a ≤ 0 || return false
  norm((p - b) × a) / norm(a) ≤ c.radius
end

function Base.in(p::Point{3,T}, t::Torus{T}) where {T}
  c, n⃗ = t.center, t.normal
  R, r = t.major, t.minor
  Q = rotation_between(n⃗, Vec{3,T}(0, 0, 1))
  x, y, z = Q * (p - c)
  (R - √(x^2 + y^2))^2 + z^2 ≤ r^2
end

function Base.in(p::Point{2}, t::Triangle{2})
  # given coordinates
  a, b, c = t.vertices
  x₁, y₁ = coordinates(a)
  x₂, y₂ = coordinates(b)
  x₃, y₃ = coordinates(c)
  x, y = coordinates(p)

  # barycentric coordinates
  λ₁ = ((y₂ - y₃) * (x - x₃) + (x₃ - x₂) * (y - y₃)) / ((y₂ - y₃) * (x₁ - x₃) + (x₃ - x₂) * (y₁ - y₃))
  λ₂ = ((y₃ - y₁) * (x - x₃) + (x₁ - x₃) * (y - y₃)) / ((y₂ - y₃) * (x₁ - x₃) + (x₃ - x₂) * (y₁ - y₃))
  λ₃ = 1 - λ₁ - λ₂

  # barycentric check
  0 ≤ λ₁ ≤ 1 && 0 ≤ λ₂ ≤ 1 && 0 ≤ λ₃ ≤ 1
end

function Base.in(p::Point{3}, t::Triangle{3})
  # given coordinates
  a, b, c = t.vertices

  # evaluate vectors defining geometry
  v₁ = b - a
  v₂ = c - a
  v₃ = p - a

  # calculate required dot products
  d₁₁ = v₁ ⋅ v₁
  d₁₂ = v₁ ⋅ v₂
  d₂₂ = v₂ ⋅ v₂
  d₃₁ = v₃ ⋅ v₁
  d₃₂ = v₃ ⋅ v₂

  # calculate reused denominator
  d = d₁₁ * d₂₂ - d₁₂ * d₁₂

  # barycentric coordinates
  λ₂ = (d₂₂ * d₃₁ - d₁₂ * d₃₂) / d
  λ₃ = (d₁₁ * d₃₂ - d₁₂ * d₃₁) / d

  # barycentric check
  λ₂ ≥ 0 && λ₃ ≥ 0 && (λ₂ + λ₃) ≤ 1
end

Base.in(p::Point, ngon::Ngon) = any(Δ -> p ∈ Δ, simplexify(ngon))

function Base.in(p::Point, poly::PolyArea)
  r = rings(poly)
  inside = sideof(p, first(r)) == :INSIDE
  if hasholes(poly)
    outside = all(sideof(p, r[i]) == :OUTSIDE for i in 2:length(r))
    inside && outside
  else
    inside
  end
end

Base.in(p::Point, m::Multi) = any(g -> p ∈ g, m.geoms)

"""
    point ∈ domain

Tells whether or not the `point` is in the `domain`.
"""
Base.in(p::Point, d::Domain) = any(e -> p ∈ e, d)
