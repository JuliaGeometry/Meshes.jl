# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Mesh{M,CRS,TP}

A mesh of geometries in a given manifold `M` with point coordinates specified
in a coordinate reference system `CRS`. Unlike a general domain, a mesh has a
well-defined topology `TP`.
"""
abstract type Mesh{M,CRS,TP} <: Domain{M,CRS} end

"""
    vertex(mesh, ind)

Return the vertex of a `mesh` at index `ind`.
"""
function vertex end

"""
    vertices(mesh)

Return the vertices of the `mesh`.
"""
vertices(m::Mesh) = [vertex(m, ind) for ind in 1:nvertices(m)]

"""
    nvertices(mesh)

Return the number of vertices of the `mesh`.
"""
nvertices(m::Mesh) = nvertices(topology(m))

"""
    faces(mesh, rank)

Return an iterator with `rank`-faces of the `mesh`.

## Examples

Consider a mesh of tetrahedrons embedded in a 3D space. We can loop over
all 3-faces (a.k.a. elements) or over all 2-faces to handle the interfaces
(i.e. triangles) between adjacent elements:

```julia
tetrahedrons = faces(mesh, 3)
triangles    = faces(mesh, 2)
segments     = faces(mesh, 1)
```
"""
faces(m::Mesh, rank) = (materialize(f, vertices(m)) for f in faces(topology(m), rank))

"""
    nfaces(mesh, rank)

Return the number of `rank`-faces of the `mesh`.
"""
nfaces(m::Mesh, rank) = nfaces(topology(m), rank)

"""
    elements(mesh)

Return the top-faces (a.k.a. elements) of the `mesh`.

## Examples

The elements of a volume embedded in 3D space can be tetrahedrons, hexahedrons,
or any 3-face. The elements of a surface embedded in 3D space can be triangles,
quadrangles or any 2-face.
"""
elements(m::Mesh) = (element(m, ind) for ind in 1:nelements(m))

"""
    element(mesh, ind)

Return the element of the `mesh` at index `ind`.
"""
element(m::Mesh, ind) = materialize(element(topology(m), ind), vertices(m))

"""
    nelements(mesh)

Return the number of elements of the `mesh`.
"""
nelements(m::Mesh) = nelements(topology(m))

"""
    facets(mesh)

Return the (top-1)-faces (a.k.a. facets) of the `mesh`.
"""
facets(m::Mesh) = (facet(m, ind) for ind in 1:nfacets(m))

"""
    facet(mesh, ind)

Return the facet of the `mesh` at index `ind`.
"""
facet(m::Mesh, ind) = materialize(facet(topology(m), ind), vertices(m))

"""
    nfacets(mesh)

Return the number of facets of the `mesh`.
"""
nfacets(m::Mesh) = nfacets(topology(m))

"""
    topoconvert(T, mesh)

Convert underlying topology of the `mesh` to topology of type `T`.

## Examples

Convert underlying topology to [`HalfEdgeTopology`](@ref) for
efficient topological relations.

```julia
newmesh = topoconvert(HalfEdgeTopology, mesh)
```
"""
topoconvert(TP::Type{<:Topology}, m::Mesh) = SimpleMesh(vertices(m), convert(TP, topology(m)))

==(m₁::Mesh, m₂::Mesh) = vertices(m₁) == vertices(m₂) && topology(m₁) == topology(m₂)

function Base.show(io::IO, ::MIME"text/plain", m::Mesh)
  t = topology(m)
  verts = vertices(m)
  elems = elements(t)
  nvert = nvertices(m)
  nelms = nelements(m)
  summary(io, m)
  println(io)
  println(io, "  $nvert vertices")
  printelms(io, verts, "  ")
  println(io)
  println(io, "  $nelms elements")
  printitr(io, elems, "  ")
end

"""
    Grid{M,CRS,Dim}

A grid of geometries in a given manifold `M` with points coordinates specified
in a coordinate reference system `CRS`, which is embedded in `Dim` dimensions.
"""
const Grid{M,CRS,Dim} = Mesh{M,CRS,GridTopology{Dim}}

"""
    vertex(grid, ijk)

Convert Cartesian index `ijk` to vertex on `grid`.
"""
vertex(g::Grid, ijk::CartesianIndex) = vertex(g, ijk.I)

"""
    xyz(grid)

Returns the coordinate vectors of each dimension of the `grid`, e.g `(x, y, z, ...)`.
The vertex `i,j,k,...` is constructed with `Point(x[i], y[j], z[k], ...)`.
"""
function xyz end

"""
    XYZ(grid)

Returns the coordinate arrays of each dimension of the `grid`, e.g `(X, Y, Z, ...)`.
The vertex `i,j,k,...` is constructed with `Point(X[i,j,k,...], Y[i,j,k,...], Z[i,j,k,...], ...)`.
"""
function XYZ end

# ----------
# FALLBACKS
# ----------

Base.size(g::Grid) = size(topology(g))

vertex(g::Grid, ind::Int) = vertex(g, CartesianIndices(size(g) .+ 1)[ind])

vertex(g::Grid, ijk::Dims) = vertex(g, LinearIndices(size(g) .+ 1)[ijk...])

Base.minimum(g::Grid) = vertex(g, ntuple(i -> 1, embeddim(g)))
Base.maximum(g::Grid) = vertex(g, size(g) .+ 1)
Base.extrema(g::Grid) = minimum(g), maximum(g)

function element(g::Grid, ind::Int)
  elem = element(topology(g), ind)
  type = pltype(elem)
  einds = indices(elem)
  cinds = CartesianIndices(size(g) .+ 1)
  verts = ntuple(i -> vertex(g, cinds[einds[i]]), nvertices(type))
  type(verts)
end

Base.eltype(g::Grid) = typeof(first(g))

Base.getindex(g::Grid, ind::Int) = element(g, ind)

Base.getindex(g::Grid, inds::AbstractVector) = [element(g, ind) for ind in inds]

Base.getindex(g::Grid, ijk::Int...) = element(g, LinearIndices(size(g))[ijk...])

@propagate_inbounds function Base.getindex(g::Grid, ijk...)
  dims = size(g)
  ranges = ntuple(i -> _asrange(dims[i], ijk[i]), embeddim(g))
  getindex(g, CartesianIndices(ranges))
end

function Base.getindex(g::Grid, I::CartesianIndices)
  @boundscheck _checkbounds(g, I)
  dims = size(I)
  odims = size(g)
  cinds = first(I):CartesianIndex(Tuple(last(I)) .+ 1)
  inds = vec(LinearIndices(odims .+ 1)[cinds])
  periodic = isperiodic(topology(g)) .&& dims .== odims
  SimpleMesh(vertices(g)[inds], GridTopology(dims, periodic))
end

_asrange(::Int, r::UnitRange{Int}) = r
_asrange(d::Int, ::Colon) = 1:d
_asrange(::Int, i::Int) = i:i

function _checkbounds(g, I)
  dims = size(g)
  ranges = I.indices
  if !all(first(r) ≥ 1 && last(r) ≤ d for (d, r) in zip(dims, ranges))
    throw(BoundsError(g, ranges))
  end
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("mesh/cartesiangrid.jl")
include("mesh/rectilineargrid.jl")
include("mesh/structuredgrid.jl")
include("mesh/simplemesh.jl")
include("mesh/transformedmesh.jl")
