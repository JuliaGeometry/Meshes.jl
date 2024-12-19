# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Polytope{K,M,CRS}

We say that a geometry is a K-polytope when it is a collection of "flat" sides
that constitute a `K`-dimensional subspace. They are called chain, polygon and
polyhedron respectively for 1D (`K=1`), 2D (`K=2`) and 3D (`K=3`) subspaces.
The parameter `K` is also known as the rank or parametric dimension
of the polytope (<https://en.wikipedia.org/wiki/Abstract_polytope>).

The term polytope expresses a particular combinatorial structure. A polyhedron,
for example, can be decomposed into faces. Each face can then be decomposed into
edges, and edges into vertices. Some conventions act as a mapping between vertices
and higher dimensional features (edges, faces, cells...), removing the need to
store all features.

Additionally, the following property must hold in order for a geometry to be considered
a polytope: the boundary of a (K+1)-polytope is a collection of K-polytopes, which may
have (K-1)-polytopes in common. See <https://en.wikipedia.org/wiki/Polytope>.

### Notes

- Type aliases are `Chain`, `Polygon`, `Polyhedron`.
"""
abstract type Polytope{K,M<:Manifold,C<:CRS} <: Geometry{M,C} end

# heper macro to define polytopes
macro polytope(type, K, N)
  structexpr = if K == 3
    quote
      struct $type{C<:CRS,Mₚ<:Manifold} <: Polytope{$K,𝔼{3},C}
        vertices::SVector{$N,Point{Mₚ,C}}
      end
    end
  else
    quote
      struct $type{M<:Manifold,C<:CRS} <: Polytope{$K,M,C}
        vertices::SVector{$N,Point{M,C}}
      end
    end
  end

  expr = quote
    $Base.@__doc__ $structexpr

    $type(vertices::NTuple{$N,P}) where {P<:Point} = $type(SVector(vertices))
    $type(vertices::Vararg{Tuple,$N}) = $type(Point.(vertices))
    $type(vertices::Vararg{P,$N}) where {P<:Point} = $type(vertices)
  end

  esc(expr)
end

# -------------------
# 1-POLYTOPE (CHAIN)
# -------------------

"""
    Chain{M,CRS}

A chain is a 1-polytope, i.e. a polytope with parametric dimension 1.
See <https://en.wikipedia.org/wiki/Polygonal_chain>.

See also [`Segment`](@ref), [`Rope`](@ref), [`Ring`](@ref).
"""
const Chain = Polytope{1}

"""
    segments(chain)

Return the segments linking consecutive points of the `chain`.
"""
function segments(c::Chain)
  v = vertices(c)
  n = length(v) - !isclosed(c)
  @inbounds (Segment(v[i], v[i + 1]) for i in 1:n)
end

"""
    close(chain)

Close the `chain`, i.e. add a segment going from the last to the first vertex.
"""
function Base.close(::Chain) end

"""
    open(chain)

Open the `chain`, i.e. remove the segment going from the last to the first vertex.
"""
function Base.open(::Chain) end

"""
    unique!(chain)

Remove duplicate vertices in the `chain`.
Closed chains remain closed.
"""
function Base.unique!(c::Chain)
  # sort vertices lexicographically
  verts = vertices(open(c))
  perms = sortperm(to.(verts))

  # remove true duplicates
  keep = Int[]
  sorted = @view verts[perms]
  for i in 1:(length(sorted) - 1)
    if !isapprox(sorted[i], sorted[i + 1])
      # save index in the original vector
      push!(keep, perms[i])
    end
  end
  push!(keep, last(perms))

  # preserve chain order
  sort!(keep)

  # update vertices in place
  copy!(verts, verts[keep])

  c
end

"""
    reverse!(chain)

Reverse the `chain` vertices in place.
"""
function Base.reverse!(::Chain) end

"""
    reverse(chain)

Reverse the `chain` vertices.
"""
Base.reverse(c::Chain) = reverse!(deepcopy(c))

"""
    angles(chain)

Return angles `∠(vᵢ-₁, vᵢ, vᵢ+₁)` at all vertices
`vᵢ` of the `chain`. If the chain is open, the first
and last vertices have no angles. Positive angles
represent a CCW rotation whereas negative angles
represent a CW rotation. In either case, the
absolute value of the angles returned is never
greater than `π`.
"""
function angles(c::Chain)
  vs = vertices(c)
  i1 = firstindex(vs) + !isclosed(c)
  i2 = lastindex(vs) - !isclosed(c)
  map(i -> ∠(vs[i - 1], vs[i], vs[i + 1]), i1:i2)
end

function (c::Chain)(t)
  if t < 0 || t > 1
    throw(DomainError(t, "c(t) is not defined for t outside [0, 1]."))
  end
  segs = segments(c)
  sums = cumsum(measure.(segs))
  sums /= last(sums)
  # find k such that sums[k] ≤ t < sums[k + 1]
  k = searchsortedfirst(sums, t) - 1
  # select segment s at index k
  s, _ = iterate(segs, k)
  # reparametrization of t within s
  ∑ₖ = iszero(k) ? zero(eltype(sums)) : sums[k]
  ∑ₖ₊₁ = sums[k + 1]
  s((t - ∑ₖ) / (∑ₖ₊₁ - ∑ₖ))
end

# implementations of Chain
include("polytopes/segment.jl")
include("polytopes/rope.jl")
include("polytopes/ring.jl")

# ---------------------
# 2-POLYTOPE (POLYGON)
# ---------------------

"""
    Polygon{M,CRS}

A polygon is a 2-polytope, i.e. a polytope with parametric dimension 2.

See also [`Ngon`](@ref) and [`PolyArea`](@ref).
"""
const Polygon = Polytope{2}

"""
    ≗(polygon₁, polygon₂)

Tells whether or not the `polygon₁` and `polygon₂`
are equal regardless of circular shifts.
"""
function ≗(p₁::Polygon, p₂::Polygon)
  rings₁ = rings(p₁)
  rings₂ = rings(p₂)
  nring₁ = length(rings₁)
  nring₂ = length(rings₂)
  nring₁ == nring₂ || return false
  all(r₁ ≗ r₂ for (r₁, r₂) in zip(rings₁, rings₂))
end

"""
    rings(polygon)

Return the outer and inner rings of the polygon.
"""
function rings end

# implementations of Polygon
include("polytopes/ngon.jl")
include("polytopes/polyarea.jl")

# ------------------------
# 3-POLYTOPE (POLYHEDRON)
# ------------------------

"""
    Polyhedron{M,CRS}

A polyhedron is a 3-polytope, i.e. a polytope with parametric dimension 3.

See also [`Tetrahedron`](@ref), [`Hexahedron`](@ref) and [`Pyramid`](@ref).
"""
const Polyhedron = Polytope{3}

# implementations of Polyhedron
include("polytopes/tetrahedron.jl")
include("polytopes/hexahedron.jl")
include("polytopes/pyramid.jl")
include("polytopes/wedge.jl")

# -----------------------
# N-POLYTOPE (FALLBACKS)
# -----------------------

"""
    paramdim(polytope)

Return the parametric dimension or rank of the polytope.
"""
paramdim(::Type{<:Polytope{K}}) where {K} = K

"""
    vertex(polytope, ind)

Return the vertex of a `polytope` at index `ind`.
"""
vertex(p::Polytope, ind) = vertices(p)[ind]

"""
    vertices(polytope)

Return the vertices of a `polytope`.
"""
vertices(p::Polytope) = p.vertices

"""
    nvertices(polytope)

Return the number of vertices of the `polytope`.
"""
nvertices(p::Polytope) = nvertices(typeof(p))

"""
    eachvertex(polytope)

Return an iterator for the vertices of the `polytope`.
"""
eachvertex(p::Polytope) = (vertex(p, i) for i in 1:nvertices(p))

"""
    unique(polytope)

Return a new `polytope` without duplicate vertices.
"""
Base.unique(p::Polytope) = unique!(deepcopy(p))

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, p::Polytope)
  name = prettyname(p)
  print(io, "$name(")
  printverts(io, vertices(p))
  print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", p::Polytope)
  summary(io, p)
  println(io)
  printelms(io, vertices(p))
end
