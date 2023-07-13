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
struct Ngon{N,Dim,T} <: Polygon{Dim,T}
  vertices::NTuple{N,Point{Dim,T}}
end

# function Ngon{N,Dim,T}(vertices::AbstractVector{Point{Dim,T}}) where {N,Dim,T}
#   P = length(vertices)
#   P == N || throw(ArgumentError("Invalid number of vertices for Ngon. Expected $N, got $P."))
#   v = ntuple(i -> @inbounds(vertices[i]), N)
#   Ngon{N,Dim,T}(v)
# end

Ngon(vertices::Vararg{Tuple,N}) where {N} = Ngon(Point.(vertices))
Ngon(vertices::Vararg{Point{Dim,T},N}) where {N,Dim,T} = Ngon{N,Dim,T}(vertices)

Ngon{N}(vertices::Vararg{Tuple,N}) where {N} = Ngon(Point.(vertices))
Ngon{N}(vertices::Vararg{Point{Dim,T},N}) where {N,Dim,T} = Ngon{N,Dim,T}(vertices)
Ngon{N}(vertices::NTuple{N,Point{Dim,T}}) where {N,Dim,T} = Ngon{N,Dim,T}(vertices)

# Ngon(vertices::AbstractVector{<:Tuple}) = Ngon(Point.(vertices))
# Ngon(vertices::AbstractVector{Point{Dim,T}}) where {Dim,T} = Ngon{length(vertices)}(vertices)
# Ngon{N}(vertices::AbstractVector{<:Tuple}) where {N} = Ngon{N}(Point.(vertices))
# Ngon{N}(vertices::AbstractVector{Point{Dim,T}}) where {N,Dim,T} = Ngon{N,Dim,T}(vertices)

# type aliases for convenience
const Triangle = Ngon{3}
const Quadrangle = Ngon{4}
const Pentagon = Ngon{5}
const Hexagon = Ngon{6}
const Heptagon = Ngon{7}
const Octagon = Ngon{8}
const Nonagon = Ngon{9}
const Decagon = Ngon{10}

issimple(::Type{<:Ngon}) = true

hasholes(::Ngon) = false

Base.unique!(ngon::Ngon) = ngon

nvertices(::Type{<:Ngon{N}}) where {N} = N

rings(ngon::Ngon) = [Ring(pointify(ngon))]

angles(ngon::Ngon) = angles(boundary(ngon))

innerangles(ngon::Ngon) = innerangles(boundary(ngon))

signarea(ngon::Ngon) = sum(signarea, simplexify(ngon))

Base.in(p::Point, ngon::Ngon) = any(Δ -> p ∈ Δ, simplexify(ngon))

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Ngon{N,Dim,T}}) where {N,Dim,T}
  v = ntuple(i -> rand(rng, Point{Dim,T}), N)
  Ngon(v)
end

# ----------
# TRIANGLES
# ----------

issimplex(::Type{<:Triangle}) = true

isconvex(::Type{<:Triangle}) = true

# specialize for performance
isconvex(::Triangle) = true

isparametrized(::Type{<:Triangle}) = true

function signarea(t::Triangle{2})
  v = t.vertices
  signarea(v[1], v[2], v[3])
end

measure(t::Triangle{2}) = abs(signarea(t))

function measure(t::Triangle{3})
  A, B, C = t.vertices
  norm((B - A) × (C - A)) / 2
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

function normal(t::Triangle{3})
  a, b, c = t.vertices
  n = (b - a) × (c - a)
  n / norm(n)
end

function (t::Triangle)(u, v)
  w = (1 - u - v)
  if (u < 0 || u > 1) || (v < 0 || v > 1) || (w < 0 || w > 1)
    throw(DomainError((u, v), "invalid barycentric coordinates for triangle."))
  end
  v₁, v₂, v₃ = coordinates.(t.vertices)
  Point(v₁ * w + v₂ * u + v₃ * v)
end

# ------------
# QUADRANGLES
# ------------

isperiodic(::Type{<:Quadrangle}) = (false, false)

isparametrized(::Type{<:Quadrangle}) = true

# Coons patch https://en.wikipedia.org/wiki/Coons_patch
function (q::Quadrangle)(u, v)
  if (u < 0 || u > 1) || (v < 0 || v > 1)
    throw(DomainError((u, v), "q(u, v) is not defined for u, v outside [0, 1]²."))
  end
  c₀₀, c₀₁, c₁₁, c₁₀ = coordinates.(q.vertices)
  Point(c₀₀ * (1 - u) * (1 - v) + c₀₁ * u * (1 - v) + c₁₀ * (1 - u) * v + c₁₁ * u * v)
end
