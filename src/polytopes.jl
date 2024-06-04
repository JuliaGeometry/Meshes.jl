# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Polytope{K,Dim,CRS}

We say that a geometry is a K-polytope when it is a collection of "flat" sides
that constitute a `K`-dimensional subspace. They are called chain, polygon and
polyhedron respectively for 1D (`K=1`), 2D (`K=2`) and 3D (`K=3`) subspaces,
embedded in a `Dim`-dimensional space with given coordinate reference system `CRS`. 
The parameter `K` is also known as the rank or parametric dimension 
of the polytope: <https://en.wikipedia.org/wiki/Abstract_polytope>.

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
abstract type Polytope{K,Dim,CRS} <: Geometry{Dim,CRS} end

# heper macro to define polytopes
macro polytope(type, K, N)
  expr = quote
    $Base.@__doc__ struct $type{Dim,C<:CRS} <: Polytope{$K,Dim,C}
      vertices::NTuple{$N,Point{Dim,C}}
    end

    $type(vertices::Vararg{Tuple,$N}) = $type(Point.(vertices))
    $type(vertices::Vararg{P,$N}) where {P<:Point} = $type(vertices)
  end
  esc(expr)
end

# -------------------
# 1-POLYTOPE (CHAIN)
# -------------------

"""
    Chain{Dim,CRS}

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

# implementations of Chain
include("polytopes/segment.jl")
include("polytopes/rope.jl")
include("polytopes/ring.jl")

# ---------------------
# 2-POLYTOPE (POLYGON)
# ---------------------

"""
    Polygon{Dim,CRS}

A polygon is a 2-polytope, i.e. a polytope with parametric dimension 2.

See also [`Ngon`](@ref) and [`PolyArea`](@ref).
"""
const Polygon = Polytope{2}

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
    Polyhedron{Dim,CRS}

A polyhedron is a 3-polytope, i.e. a polytope with parametric dimension 3.

See also [`Tetrahedron`](@ref), [`Hexahedron`](@ref) and [`Pyramid`](@ref).
"""
const Polyhedron = Polytope{3}

# implementations of Polyhedron
include("polytopes/tetrahedron.jl")
include("polytopes/hexahedron.jl")
include("polytopes/pyramid.jl")

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

Return the number of vertices in the `polytope`.
"""
nvertices(p::Polytope) = nvertices(typeof(p))

"""
    centroid(polytope)

Return the centroid of the `polytope`.
"""
centroid(p::Polytope) = Point(coords(sum(to, vertices(p)) / length(vertices(p))))

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
