# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    point ∈ geometry

Tells whether or not the `point` is in the `geometry`.
"""
Base.in(p::Point, g::Geometry) = sideof(p, boundary(g)) == IN

Base.in(p₁::Point, p₂::Point) = p₁ == p₂

function Base.in(p::Point, s::Segment)
  # given collinear points (a, b, p), the point p intersects
  # segment ab if and only if vectors satisfy 0 ≤ ap ⋅ ab ≤ ||ab||²
  a, b = vertices(s)
  ab, ap = b - a, p - a
  iscollinear(a, b, p) && (abap = ab ⋅ ap;
  isnonnegative(abap) && abap ≤ ab ⋅ ab)
end

Base.in(p::Point, r::Ray) = p ∈ Line(r(0), r(1)) && isnonnegative((p - r(0)) ⋅ (r(1) - r(0)))

function Base.in(p::Point, l::Line)
  w = norm(l(1) - l(0))
  d = evaluate(Euclidean(), p, l)
  d + w ≈ w # d ≈ 0.0 will be too precise, and d < atol{T} can't scale.
end

Base.in(p::Point, c::Chain) = any(s -> p ∈ s, segments(c))

Base.in(p::Point, pl::Plane) = isapproxzero(udot(normal(pl), p - pl(0, 0)))

Base.in(p::Point, b::Box) = minimum(b) ⪯ p ⪯ maximum(b)

function Base.in(p::Point, b::Ball)
  c = center(b)
  r = radius(b)
  s = norm(p - c)
  s < r || isapproxequal(s, r)
end

function Base.in(p::Point, s::Sphere)
  c = center(s)
  r = radius(s)
  s = norm(p - c)
  isapproxequal(s, r)
end

function Base.in(p::Point, d::Disk)
  p ∉ plane(d) && return false
  c = center(d)
  r = radius(d)
  s = norm(p - c)
  s < r || isapproxequal(s, r)
end

function Base.in(p::Point, c::Circle)
  p ∉ plane(c) && return false
  o = center(c)
  r = radius(c)
  s = norm(p - o)
  isapproxequal(s, r)
end

function Base.in(p::Point, c::Cone)
  a = apex(c)
  b = center(base(c))
  ax = a - b
  isnonnegative((a - p) ⋅ ax) || return false
  isnonpositive((b - p) ⋅ ax) || return false
  ∠(b, a, p) ≤ halfangle(c)
end

function Base.in(p::Point, c::Cylinder)
  b = bottom(c)(0, 0)
  t = top(c)(0, 0)
  r = radius(c)
  a = t - b
  isnonnegative((p - b) ⋅ a) || return false
  isnonpositive((p - t) ⋅ a) || return false
  norm((p - b) × a) / norm(a) ≤ r
end

function Base.in(p::Point, f::Frustum)
  t = center(top(f))
  b = center(bottom(f))
  ax = b - t
  isnonnegative((p - t) ⋅ ax) || return false
  isnonpositive((p - b) ⋅ ax) || return false
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

function Base.in(p::Point, t::Torus)
  ℒ = lentype(p)
  R, r = radii(t)
  c, n = center(t), normal(t)
  Q = urotbetween(n, Vec(zero(ℒ), zero(ℒ), oneunit(ℒ)))
  x, y, z = Q * (p - c)
  (R - √(x^2 + y^2))^2 + z^2 ≤ r^2
end

function Base.in(point::Point, poly::Polygon{𝔼{2}})
  r = rings(poly)
  inside = sideof(point, first(r)) != OUT
  if hasholes(poly)
    outside = all(sideof(point, r[i]) == OUT for i in 2:length(r))
    inside && outside
  else
    inside
  end
end

Base.in(p::Point, poly::Polygon{𝔼{3}}) = any(Δ -> p ∈ Δ, simplexify(poly))

function Base.in(p::Point, t::Triangle{𝔼{3}})
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

Base.in(p::Point, m::Multi) = any(g -> p ∈ g, parent(m))

"""
    point ∈ domain

Tells whether or not the `point` is in the `domain`.
"""
Base.in(p::Point, d::Domain) = any(e -> p ∈ e, d)
