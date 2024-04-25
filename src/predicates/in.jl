# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    point ∈ geometry

Tells whether or not the `point` is in the `geometry`.
"""
Base.in(p::Point, g::Geometry) = sideof(p, boundary(g)) == IN

Base.in(p₁::Point, p₂::Point) = p₁ == p₂

function Base.in(p::Point{Dim}, s::Segment{Dim}) where {Dim}
  # given collinear points (a, b, p), the point p intersects
  # segment ab if and only if vectors satisfy 0 ≤ ap ⋅ ab ≤ ||ab||²
  a, b = vertices(s)
  ab, ap = b - a, p - a
  iscollinear(a, b, p) && zero(coordtype(p)) ≤ ab ⋅ ap ≤ ab ⋅ ab
end

Base.in(p::Point, r::Ray) = p ∈ Line(r(0), r(1)) && (p - r(0)) ⋅ (r(1) - r(0)) ≥ 0

function Base.in(p::Point, l::Line)
  w = norm(l(1) - l(0))
  d = evaluate(Euclidean(), p, l)
  d + w ≈ w # d ≈ 0.0 will be too precise, and d < atol{T} can't scale.
end

Base.in(p::Point, c::Chain) = any(s -> p ∈ s, segments(c))

function Base.in(p::Point{3}, pl::Plane)
  𝒬 = coordtype(p)
  isapprox(normal(pl) ⋅ (p - pl(0, 0)), zero(𝒬), atol=atol(𝒬))
end

Base.in(p::Point, b::Box) = minimum(b) ⪯ p ⪯ maximum(b)

function Base.in(p::Point{Dim}, b::Ball{Dim}) where {Dim}
  c = center(b)
  r = radius(b)
  s = norm(p - c)
  s < r || isapprox(s, r, atol=atol(coordtype(p)))
end

function Base.in(p::Point{Dim}, s::Sphere{Dim}) where {Dim}
  c = center(s)
  r = radius(s)
  s = norm(p - c)
  isapprox(s, r, atol=atol(coordtype(p)))
end

function Base.in(p::Point{3}, d::Disk)
  p ∉ plane(d) && return false
  c = center(d)
  r = radius(d)
  s = norm(p - c)
  s < r || isapprox(s, r, atol=atol(coordtype(p)))
end

function Base.in(p::Point{3}, c::Circle)
  p ∉ plane(c) && return false
  o = center(c)
  r = radius(c)
  s = norm(p - o)
  isapprox(s, r, atol=atol(coordtype(p)))
end

function Base.in(p::Point{3}, c::Cone)
  a = apex(c)
  b = center(base(c))
  ax = a - b
  (a - p) ⋅ ax ≥ 0 || return false
  (b - p) ⋅ ax ≤ 0 || return false
  ∠(b, a, p) ≤ halfangle(c)
end

function Base.in(p::Point{3}, c::Cylinder)
  b = bottom(c)(0, 0)
  t = top(c)(0, 0)
  r = radius(c)
  a = t - b
  (p - b) ⋅ a ≥ 0 || return false
  (p - t) ⋅ a ≤ 0 || return false
  norm((p - b) × a) / norm(a) ≤ r
end

function Base.in(p::Point{3}, f::Frustum)
  t = center(top(f))
  b = center(bottom(f))
  ax = b - t
  (p - t) ⋅ ax ≥ 0 || return false
  (p - b) ⋅ ax ≤ 0 || return false
  # axial distance of p
  ad = (p - t) ⋅ normalize(ax)
  adrel = ad / norm(ax)
  # frustum radius at axial distance of p
  rt = radius(top(f))
  rb = radius(bottom(f))
  r = rt * (1 - adrel) + rb * adrel
  # radial distance of p
  rd = norm((p - t) - adrel * ax)
  rd ≤ r
end

function Base.in(p::Point{3}, t::Torus)
  𝒬 = coordtype(p)
  R, r = radii(t)
  c, n = center(t), normal(t)
  Q = rotation_between(n, Vec(zero(𝒬), zero(𝒬), one(𝒬)))
  x, y, z = Q * (p - c)
  (R - √(x^2 + y^2))^2 + z^2 ≤ r^2
end

function Base.in(p::Point{2}, t::Triangle{2})
  # given coordinates
  a, b, c = vertices(t)
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
  a, b, c = vertices(t)

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
  inside = sideof(p, first(r)) == IN
  if hasholes(poly)
    outside = all(sideof(p, r[i]) == OUT for i in 2:length(r))
    inside && outside
  else
    inside
  end
end

Base.in(p::Point, m::Multi) = any(g -> p ∈ g, parent(m))

"""
    point ∈ domain

Tells whether or not the `point` is in the `domain`.
"""
Base.in(p::Point, d::Domain) = any(e -> p ∈ e, d)
