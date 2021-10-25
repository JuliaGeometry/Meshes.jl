# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Polytope{K,Dim,T}

We say that a geometry is a K-polytope when it is a collection of "flat" sides
that constitue a `K`-dimensional subspace. They are called polygon and polyhedron
respectively for 2D (`K=2`) and 3D (`K=3`) subspaces, embedded in a `Dim`-dimensional
space. The parameter `K` is also known as the rank or parametric dimension of the
polytope: https://en.wikipedia.org/wiki/Abstract_polytope.

The term polytope expresses a particular combinatorial structure. A polyhedron,
for example, can be decomposed into faces. Each face can then be decomposed into
edges, and edges into vertices. Some conventions act as a mapping between vertices
and higher dimensional features (edges, faces, cells...), removing the need to
store all features.

Additionally, the following property must hold in order for a geometry to be considered
a polytope: the boundary of a (K+1)-polytope is a collection of K-polytopes, which may
have (K-1)-polytopes in common. See https://en.wikipedia.org/wiki/Polytope.

### Notes

- Type aliases are `Polygon`, `Polyhedron`.
"""
abstract type Polytope{K,Dim,T} <: Geometry{Dim,T} end

(::Type{PL})(vertices::Vararg{P}) where {PL<:Polytope,P<:Point} = PL(SVector(vertices))
(::Type{PL})(vertices::AbstractVector{TP}) where {PL<:Polytope,TP<:Tuple} = PL(Point.(vertices))
(::Type{PL})(vertices::Vararg{TP}) where {PL<:Polytope,TP<:Tuple} = PL(collect(vertices))

"""
    paramdim(polytope)

Return the parametric dimension or rank of the polytope.
"""
paramdim(::Type{<:Polytope{K}}) where {K} = K

"""
    vertices(polytope)

Return the vertices of a `polytope`.
"""
vertices(p::Polytope) = p.vertices

"""
    nvertices(polytope)

Return the number of vertices in the `polytope`.
"""
nvertices(p::Polytope) = length(vertices(p))

"""
    p1 == p2

Tells whether or not polytopes `p1` and `p2` are equal.
"""
==(p1::Polytope, p2::Polytope) = vertices(p1) == vertices(p2)

"""
    centroid(polytope)

Return the centroid of the `polytope`.
"""
centroid(p::Polytope) = Point(sum(coordinates, vertices(p)) / length(vertices(p)))

"""
    unique(polytope)

Return a new `polytope` without duplicate vertices.
"""
Base.unique(p::Polytope) = unique!(deepcopy(p))

"""
    unique!(polytope)

Remove duplicate vertices in `polytope`.
"""
function Base.unique!(::Polytope) end

function Base.show(io::IO, p::Polytope)
  kind = prettyname(p)
  vert = join(vertices(p), ", ")
  print(io, "$kind($vert)")
end

function Base.show(io::IO, ::MIME"text/plain", p::Polytope{K,Dim,T}) where {K,Dim,T}
  kind = prettyname(p)
  lines = ["  └─$v" for v in vertices(p)]
  println(io, "$kind{$Dim,$T}")
  print(io, join(lines, "\n"))
end

prettyname(p::Polytope) = prettyname(typeof(p))
function prettyname(PL::Type{<:Polytope})
  n = string(PL)
  i = findfirst('{', n)
  isnothing(i) ? n : n[1:i-1]
end

# --------
# POLYGON
# --------

"""
    Polygon{Dim,T}

A polygon is a 2-polytope, i.e. a polytope with parametric dimension 2.
"""
const Polygon = Polytope{2}

"""
    area(polygon)

Return the area of the `polygon`.
"""
area(p::Polygon) = measure(p)

"""
    chains(polygon)

Return the outer and inner chains of the polygon.
"""
function chains(::Polygon) end

"""
    hasholes(polygon)

Tells whether or not the `polygon` contains holes.
"""
function hasholes(::Polygon) end

"""
    issimple(polygon)

Tells whether or not the `polygon` is simple.
See [https://en.wikipedia.org/wiki/Simple_polygon]
(https://en.wikipedia.org/wiki/Simple_polygon).
"""
issimple(p::Polygon) = issimple(typeof(p))

"""
    windingnumber(point, polygon)

Winding number of `point` with respect to the `polygon`.
"""
function windingnumber(::Point, ::Polygon) end

"""
    orientation(polygon)

Returns the orientation of the rings of the `polygon`
as either counter-clockwise (CCW) or clockwise (CW).
"""
orientation(p::Polygon) = orientation(p, WindingOrientation())

function orientation(p::Polygon, algo)
  o = [orientation(c, algo) for c in chains(p)]
  hasholes(p) ? o : first(o)
end

"""
    boundary(polygon)

Returns the boundary of the `polygon`.
"""
boundary(p::Polygon) = hasholes(p) ? Multi(chains(p)) : first(chains(p))

"""
    isconvex(polygon)

Tells whether or not the `polygon` is convex.
"""
isconvex(p::Polygon{Dim,T}) where {Dim,T} =
  issimple(p) && all(≤(T(π)), innerangles(boundary(p)))

"""
    bridge(polygon; width=0)

Transform `polygon` with holes into a single outer chain
via bridges of given `width` as described in Held 1998.
Return the outer chain and a vector with pairs of indices
for duplicate vertices. These indices can be used to undo
the bridges.

## References

* Held. 1998. [FIST: Fast Industrial-Strength Triangulation of Polygons]
  (https://link.springer.com/article/10.1007/s00453-001-0028-4)
"""
function bridge(p::Polygon{Dim,T}; width=zero(T)) where {Dim,T}
  # polygons without holes are trivial
  if !hasholes(p)
    outerchain = first(chains(p))
    duplicates = Tuple{Int,Int}[]
    return outerchain, duplicates
  end

  # retrieve chains as vectors of coordinates
  pchains = [coordinates.(vertices(open(c))) for c in chains(p)]

  # sort vertices lexicographically
  coords  = [coord for pchain in pchains for coord in pchain]
  indices = sortperm(sortperm(coords))

  # each chain has its own set of indices
  pinds = Vector{Int}[]; offset = 0
  for nvertex in length.(pchains)
    push!(pinds, indices[offset+1:offset+nvertex])
    offset += nvertex
  end

  # sort chains based on leftmost vertex
  leftmost = argmin.(pinds)
  minimums = getindex.(pinds, leftmost)
  reorder  = sortperm(minimums)
  leftmost = leftmost[reorder]
  minimums = minimums[reorder]
  pchains  = pchains[reorder]
  pinds    = pinds[reorder]

  # initialize outer boundary
  outer = first(pchains)
  oinds = first(pinds)

  # merge holes into outer boundary
  for i in 2:length(pchains)
    inner = pchains[i]
    iinds = pinds[i]
    l = leftmost[i]
    m = minimums[i]

    # find closest vertex in boundary
    dmin, jmin = Inf, 0
    for j in findall(oinds .≤ m)
      d = sum(abs, outer[j] - inner[l])
      if d < dmin
        dmin, jmin = d, j
      end
    end

    # create a bridge of given width δ
    # from line segment A--B. The point
    # A is split into A′ and A′′ and the
    # point B is split into B′ and B′′
    A = outer[jmin]
    B = inner[l]
    δ = width
    v = B - A
    u = Vec(-v[2], v[1])
    n = u / norm(u)
    A′  = A + δ/2 * n
    A′′ = A - δ/2 * n
    B′  = B + δ/2 * n
    B′′ = B - δ/2 * n

    # insert hole at closest vertex
    outer = [
      outer[begin:jmin-1]; [A′, B′];
      circshift(inner, -l+1)[2:end];
      [B′′, A′′]; outer[jmin+1:end]
    ]
    oinds = [
      oinds[begin:jmin];
      circshift(iinds, -l+1);
      [iinds[l]];
      oinds[jmin:end]
    ]
  end

  # find duplicate vertices
  duplicates = Tuple{Int,Int}[]
  occurred = Dict{Int,Int}()
  for (i, ind) in enumerate(oinds)
    if haskey(occurred, ind)
      push!(duplicates, (occurred[ind], i))
    else
      occurred[ind] = i
    end
  end

  # close outer boundary
  push!(outer, first(outer))

  outerchain = Chain(Point.(outer))

  outerchain, duplicates
end

# -----------
# POLYHEDRON
# -----------

"""
    Polyhedron{Dim,T}

A polyhedron is a 3-polytope, i.e. a polytope with parametric dimension 3.
"""
const Polyhedron = Polytope{3}

"""
   volume(polyhedron)

Return the volume of the `polyhedron`.
"""
volume(p::Polyhedron) = measure(p)

# ----------------
# IMPLEMENTATIONS
# ----------------

include("polytopes/segment.jl")
include("polytopes/ngon.jl")
include("polytopes/chain.jl")
include("polytopes/polyarea.jl")
include("polytopes/pyramid.jl")
include("polytopes/tetrahedron.jl")
include("polytopes/hexahedron.jl")
