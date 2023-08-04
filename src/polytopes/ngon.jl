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

Ngon(vertices::Vararg{Tuple,N}) where {N} = Ngon(Point.(vertices))
Ngon(vertices::Vararg{Point{Dim,T},N}) where {N,Dim,T} = Ngon{N,Dim,T}(vertices)

Ngon{N}(vertices::Vararg{Tuple,N}) where {N} = Ngon(Point.(vertices))
Ngon{N}(vertices::Vararg{Point{Dim,T},N}) where {N,Dim,T} = Ngon{N,Dim,T}(vertices)
Ngon{N}(vertices::NTuple{N,Point{Dim,T}}) where {N,Dim,T} = Ngon{N,Dim,T}(vertices)

# type aliases for convenience
const Triangle = Ngon{3}
const Quadrangle = Ngon{4}
const Pentagon = Ngon{5}
const Hexagon = Ngon{6}
const Heptagon = Ngon{7}
const Octagon = Ngon{8}
const Nonagon = Ngon{9}
const Decagon = Ngon{10}

Base.unique!(ngon::Ngon) = ngon

nvertices(::Type{<:Ngon{N}}) where {N} = N

rings(ngon::Ngon) = [Ring(pointify(ngon))]

angles(ngon::Ngon) = angles(boundary(ngon))

innerangles(ngon::Ngon) = innerangles(boundary(ngon))

signarea(ngon::Ngon) = sum(signarea, simplexify(ngon))

# ----------
# TRIANGLES
# ----------

function signarea(t::Triangle{2})
  v = t.vertices
  signarea(v[1], v[2], v[3])
end

measure(t::Triangle{2}) = abs(signarea(t))

function measure(t::Triangle{3})
  A, B, C = t.vertices
  norm((B - A) × (C - A)) / 2
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

# Coons patch https://en.wikipedia.org/wiki/Coons_patch
function (q::Quadrangle)(u, v)
  if (u < 0 || u > 1) || (v < 0 || v > 1)
    throw(DomainError((u, v), "q(u, v) is not defined for u, v outside [0, 1]²."))
  end
  c₀₀, c₀₁, c₁₁, c₁₀ = coordinates.(q.vertices)
  Point(c₀₀ * (1 - u) * (1 - v) + c₀₁ * u * (1 - v) + c₁₀ * (1 - u) * v + c₁₁ * u * v)
end
