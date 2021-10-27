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

nvertices(::Type{<:Ngon{N}}) where {N} = N
nvertices(ngon::Ngon) = nvertices(typeof(ngon))

# measure of N-gon embedded in 2D
function signarea(ngon::Ngon{N,2}) where {N}
  v = ngon.vertices
  sum(i -> signarea(v[1], v[i], v[i+1]), 2:N-1)
end
measure(ngon::Ngon{N,2}) where {N} = abs(signarea(ngon))

# measure of N-gon embedded in higher dimension
function measure(ngon::Ngon{N}) where {N}
  areaₜ(A, B, C) = norm((B - A) × (C - A)) / 2
  v = ngon.vertices
  sum(i -> areaₜ(v[1], v[i], v[i+1]), 2:N-1)
end

hasholes(::Ngon) = false

chains(ngon::Ngon{N}) where {N} = [Chain(ngon.vertices[[1:N; 1]])]

Base.unique!(ngon::Ngon) = ngon

"""
    angles(ngon)

Return the angles of the boundary of the `ngon`.

See also [`Chain`](@ref).
"""
angles(ngon::Ngon) = angles(boundary(ngon))

"""
    innerangles(ngon)

Return inner angles of the boundary of the `ngon`.

See also [`Chain`](@ref).
"""
innerangles(ngon::Ngon) = innerangles(boundary(ngon))
    
function Base.in(p::Point{Dim,T}, ngon::Ngon{N,Dim,T}) where {N,Dim,T}
  # decompose n-gons into triangles by
  # fan triangulation (assumes convexity)
  # https://en.wikipedia.org/wiki/Fan_triangulation
  v = ngon.vertices
  Δ(i) = Triangle(view(v, [1,i,i+1]))
  any(i -> p ∈ Δ(i), 2:N-1)
end

# ----------
# TRIANGLES
# ----------

# triangles are special
issimplex(::Type{<:Triangle}) = true
isconvex(::Type{<:Triangle}) = true

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

"""
    normal(triangle)

Determine the normalised normal of `triangle` in three
dimensions
"""
function normal(t::Triangle{3})
  a, b, c = t.vertices
  n = (b - a) × (c - a)
  n / norm(n)
end

# ------------
# QUADRANGLES
# ------------

# Coons patch https://en.wikipedia.org/wiki/Coons_patch 
function (q::Quadrangle)(u, v)
  c₀₀, c₀₁, c₁₁, c₁₀ = coordinates.(q.vertices)
  Point(c₀₀*(1-u)*(1-v) + c₀₁*u*(1-v) + c₁₀*(1-u)*v + c₁₁*u*v)
end
