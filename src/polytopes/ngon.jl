# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ngon(p₁, p₂, ..., pₙ)

A N-gon is a polygon with `N ≥ 3` vertices `p₁`, `p₂`, ..., `pₙ`
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
struct Ngon{N,Dim,C<:CRS} <: Polygon{Dim,C}
  vertices::NTuple{N,Point{Dim,C}}
  function Ngon{N,Dim,C}(vertices) where {N,Dim,C<:CRS}
    if N < 3
      throw(ArgumentError("the number of vertices must be greater than or equal to 3"))
    end
    new(vertices)
  end
end

Ngon{N}(vertices::NTuple{N,Point{Dim,C}}) where {N,Dim,C<:CRS} = Ngon{N,Dim,C}(vertices)
Ngon{N}(vertices::Vararg{P,N}) where {N,P<:Point} = Ngon{N}(vertices)
Ngon{N}(vertices::Vararg{Tuple,N}) where {N} = Ngon{N}(Point.(vertices))

Ngon(vertices::NTuple{N,Point{Dim,C}}) where {N,Dim,C<:CRS} = Ngon{N,Dim,C}(vertices)
Ngon(vertices::P...) where {P<:Point} = Ngon(vertices)
Ngon(vertices::Tuple...) = Ngon(Point.(vertices))

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

function Base.isapprox(p₁::Ngon, p₂::Ngon; kwargs...)
  nvertices(p₁) ≠ nvertices(p₂) && return false
  all(isapprox(v₁, v₂; kwargs...) for (v₁, v₂) in zip(p₁.vertices, p₂.vertices))
end

rings(ngon::Ngon) = [Ring(pointify(ngon))]

angles(ngon::Ngon) = angles(boundary(ngon))

innerangles(ngon::Ngon) = innerangles(boundary(ngon))

signarea(ngon::Ngon) = sum(signarea, simplexify(ngon))

Random.rand(rng::Random.AbstractRNG, ::Type{Ngon{N,Dim}}) where {N,Dim} = Ngon{N}(ntuple(i -> rand(rng, Point{Dim}), N))

# ----------
# TRIANGLES
# ----------

function signarea(t::Triangle{2})
  v = t.vertices
  signarea(v[1], v[2], v[3])
end

signarea(::Triangle{3}) = error("signed area only defined for triangles embedded in R², use `area` instead")

function normal(t::Triangle{3})
  A, B, C = t.vertices
  unormalize(ucross((B - A), (C - A)))
end

function (t::Triangle)(u, v)
  w = (1 - u - v)
  if (u < 0 || u > 1) || (v < 0 || v > 1) || (w < 0 || w > 1)
    throw(DomainError((u, v), "invalid barycentric coordinates for triangle."))
  end
  v₁, v₂, v₃ = to.(t.vertices)
  Point(coords(v₁ * w + v₂ * u + v₃ * v))
end

# ------------
# QUADRANGLES
# ------------

# Coons patch https://en.wikipedia.org/wiki/Coons_patch
function (q::Quadrangle)(u, v)
  if (u < 0 || u > 1) || (v < 0 || v > 1)
    throw(DomainError((u, v), "q(u, v) is not defined for u, v outside [0, 1]²."))
  end
  c₀₀, c₀₁, c₁₁, c₁₀ = to.(q.vertices)
  Point(coords(c₀₀ * (1 - u) * (1 - v) + c₀₁ * u * (1 - v) + c₁₀ * (1 - u) * v + c₁₁ * u * v))
end
