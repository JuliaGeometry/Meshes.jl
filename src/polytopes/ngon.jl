# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ngon(p1, p2, ..., pN)

A N-gon is a polygon with `N` vertices `p1`, `p2`, ..., `pN`
oriented counter-clockwise (CCW). In this case the number of
vertices is fixed and known at compile time. Examples of N-gon
are `Triangle` (N=3), `Quadrangle` (N=4), `Pentagon` (N=5), etc.

### Notes

- Although the number of vertices `N` is known at compile time,
  we use abstract vectors to store the list of vertices. This
  design allows constructing N-gon from views of global vectors
  without expensive memory allocations.

- Type aliases are `Triangle`, `Quadrangle`, `Pentagon`, `Hexagon`,
  `Heptagon`, `Octagon`, `Nonagon`, `Decagon`.
"""
struct Ngon{N,Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Polygon{Dim,T}
  vertices::V
end

Ngon{N}(vertices::AbstractVector{Point{Dim,T}}) where {N,Dim,T} =
  Ngon{N,Dim,T,typeof(vertices)}(vertices)

Ngon(vertices::AbstractVector{Point{Dim,T}}) where {Dim,T} =
  Ngon{length(vertices)}(vertices)

# type aliases for convenience
const Triangle   = Ngon{3}
const Quadrangle = Ngon{4}
const Pentagon   = Ngon{5}
const Hexagon    = Ngon{6}
const Heptagon   = Ngon{7}
const Octagon    = Ngon{8}
const Nonagon    = Ngon{9}
const Decagon    = Ngon{10}

issimple(::Type{<:Ngon}) = true

hasholes(::Ngon) = false

Base.unique!(ngon::Ngon) = ngon

nvertices(::Type{<:Ngon{N}}) where {N} = N
nvertices(ngon::Ngon) = nvertices(typeof(ngon))

chains(ngon::Ngon{N}) where {N} = [Chain(ngon.vertices[[1:N; 1]])]

angles(ngon::Ngon) = angles(boundary(ngon))

innerangles(ngon::Ngon) = innerangles(boundary(ngon))

signarea(ngon::Ngon) = sum(signarea, triangulate(ngon))

measure(ngon::Ngon) = sum(measure, triangulate(ngon))

Base.in(p::Point, ngon::Ngon) = any(Δ -> p ∈ Δ, triangulate(ngon))

# ----------
# TRIANGLES
# ----------

# triangles are special
issimplex(::Type{<:Triangle}) = true
isconvex(::Type{<:Triangle}) = true

function signarea(t::Triangle{2})
  v = t.vertices
  signarea(v[1], v[2], v[3])
end

measure(t::Triangle{2}) = abs(signarea(t))

function measure(t::Triangle{3})
  A, B, C = t.vertices
  norm((B - A) × (C - A)) / 2
end

function Base.in(p::Point{2,T}, t::Triangle{2,T}) where {T}
  # given coordinates
  a, b, c = t.vertices
  x₁, y₁ = coordinates(a)
  x₂, y₂ = coordinates(b)
  x₃, y₃ = coordinates(c)
  x , y  = coordinates(p)

  # barycentric coordinates
  λ₁ = ((y₂ - y₃)*(x  - x₃) + (x₃ - x₂)*(y  - y₃)) /
       ((y₂ - y₃)*(x₁ - x₃) + (x₃ - x₂)*(y₁ - y₃))
  λ₂ = ((y₃ - y₁)*(x  - x₃) + (x₁ - x₃)*(y  - y₃)) /
       ((y₂ - y₃)*(x₁ - x₃) + (x₃ - x₂)*(y₁ - y₃))
  λ₃ = one(T) - λ₁ - λ₂

  # barycentric check
  zero(T) ≤ λ₁ ≤ one(T) &&
  zero(T) ≤ λ₂ ≤ one(T) &&
  zero(T) ≤ λ₃ ≤ one(T)
end

function Base.in(p::Point{3,T}, t::Triangle{3,T}) where {T}
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
  (λ₂ ≥ zero(T)) && (λ₃ ≥ zero(T)) && ((λ₂ + λ₃) ≤ one(T))
end

function normal(t::Triangle{3,T}) where {T}
  a, b, c = t.vertices
  n = (b - a) × (c - a)
  n / norm(n) * oneunit(T)
end

function (t::Triangle)(u::T, v::T) where {T}
  w = (one(T) - u - v)
  if u < zero(T) || u > one(T) ||
     v < zero(T) || v > one(T) ||
     w < zero(T) || w > one(T)
     throw(DomainError("barycentric coordinates out of range"))
  end
  v₁, v₂, v₃ = coordinates.(t.vertices)
  Point(v₁*w + v₂*u + v₃*v)
end

# ------------
# QUADRANGLES
# ------------

# Coons patch https://en.wikipedia.org/wiki/Coons_patch
function (q::Quadrangle)(u, v)
  c₀₀, c₀₁, c₁₁, c₁₀ = coordinates.(q.vertices)
  Point(c₀₀*(1-u)*(1-v) + c₀₁*u*(1-v) + c₁₀*(1-u)*v + c₁₁*u*v)
end
