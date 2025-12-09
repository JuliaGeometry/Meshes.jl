# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Mesh{M,CRS,TP}

A mesh of geometries in a given manifold `M` with point coordinates specified
in a coordinate reference system `CRS`. Unlike a general domain, a mesh has a
well-defined topology `TP`.
"""
abstract type Mesh{M<:Manifold,C<:CRS,TP<:Topology} <: Domain{M,C} end

"""
    vertex(mesh, ind)

Return the vertex of a `mesh` at index `ind`.
"""
function vertex end

"""
    vertices(mesh)

Return the vertices of the `mesh`.
"""
vertices(m::Mesh) = collect(eachvertex(m))

"""
    nvertices(mesh)

Return the number of vertices of the `mesh`.
"""
nvertices(m::Mesh) = nvertices(topology(m))

"""
    eachvertex(mesh)

Return an iterator for the vertices of the `mesh`.
"""
eachvertex(m::Mesh) = (vertex(m, i) for i in 1:nvertices(m))

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
  nvert = nvertices(m)
  nelms = nelements(m)
  summary(io, m)
  println(io)
  println(io, "  $nvert vertices")
  printelms(io, m, nelms=nvert, getelm=vertex, tab="  ")
  println(io)
  println(io, "  $nelms elements")
  printelms(io, t, nelms=nelms, getelm=element, tab="  ")
end

"""
    Grid{M,CRS,N}

A `N`-dimensional grid of geometries in a given manifold `M` with
point coordinates specified in a coordinate reference system `CRS`.
"""
const Grid{M<:Manifold,C<:CRS,N} = Mesh{M,C,GridTopology{N}}

"""
    vertex(grid, ijk)

Convert Cartesian index `ijk` to vertex on `grid`.
"""
vertex(g::Grid, ijk::CartesianIndex) = vertex(g, ijk.I)

"""
    vsize(grid)

Number of vertices along each dimension of the `grid`.
"""
vsize(g::Grid) = size(g) .+ .!isperiodic(g)

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
XYZ(grid::Grid) = XYZ(grid, ntuple(i -> 1, paramdim(grid)))

XYZ(grid::Grid, factors::Dims) = XYZ(grid, Val(paramdim(grid)), factors)

function XYZ(grid::Grid, ::Val{1}, factors::Dims{1})
  T = numtype(lentype(grid))
  fᵢ = only(factors)
  sᵢ = only(size(grid))
  us = 0:T(1 / fᵢ):1
  catᵢ(A...) = cat(A..., dims=Val(1))

  xyz(p) = CoordRefSystems.values(coords(p))
  mat(seg) = [xyz(seg(u)) for u in us]
  M = [mat(grid[i]) for i in 1:sᵢ]

  C = mapreduce(catᵢ, 1:sᵢ) do i
    Mᵢ = M[i]
    i == sᵢ ? Mᵢ : Mᵢ[begin:(end - 1)]
  end

  ntuple(i -> getindex.(C, i), CoordRefSystems.ncoords(crs(grid)))
end

function XYZ(grid::Grid, ::Val{2}, factors::Dims{2})
  T = numtype(lentype(grid))
  fᵢ, fⱼ = factors
  sᵢ, sⱼ = size(grid)
  us = 0:T(1 / fᵢ):1
  vs = 0:T(1 / fⱼ):1
  catᵢ(A...) = cat(A..., dims=Val(1))
  catⱼ(A...) = cat(A..., dims=Val(2))

  xyz(p) = CoordRefSystems.values(coords(p))
  mat(quad) = [xyz(quad(u, v)) for u in us, v in vs]
  M = [mat(grid[i, j]) for i in 1:sᵢ, j in 1:sⱼ]

  C = mapreduce(catⱼ, 1:sⱼ) do j
    Mⱼ = mapreduce(catᵢ, 1:sᵢ) do i
      Mᵢⱼ = M[i, j]
      i == sᵢ ? Mᵢⱼ : Mᵢⱼ[begin:(end - 1), :]
    end
    j == sⱼ ? Mⱼ : Mⱼ[:, begin:(end - 1)]
  end

  ntuple(i -> getindex.(C, i), CoordRefSystems.ncoords(crs(grid)))
end

function XYZ(grid::Grid, ::Val{3}, factors::Dims{3})
  T = numtype(lentype(grid))
  fᵢ, fⱼ, fₖ = factors
  sᵢ, sⱼ, sₖ = size(grid)
  us = 0:T(1 / fᵢ):1
  vs = 0:T(1 / fⱼ):1
  ws = 0:T(1 / fₖ):1
  catᵢ(A...) = cat(A..., dims=Val(1))
  catⱼ(A...) = cat(A..., dims=Val(2))
  catₖ(A...) = cat(A..., dims=Val(3))

  xyz(p) = CoordRefSystems.values(coords(p))
  mat(hex) = [xyz(hex(u, v, w)) for u in us, v in vs, w in ws]
  M = [mat(grid[i, j, k]) for i in 1:sᵢ, j in 1:sⱼ, k in 1:sₖ]

  C = mapreduce(catₖ, 1:sₖ) do k
    Mₖ = mapreduce(catⱼ, 1:sⱼ) do j
      Mⱼₖ = mapreduce(catᵢ, 1:sᵢ) do i
        Mᵢⱼₖ = M[i, j, k]
        i == sᵢ ? Mᵢⱼₖ : Mᵢⱼₖ[begin:(end - 1), :, :]
      end
      j == sⱼ ? Mⱼₖ : Mⱼₖ[:, begin:(end - 1), :]
    end
    k == sₖ ? Mₖ : Mₖ[:, :, begin:(end - 1)]
  end

  ntuple(i -> getindex.(C, i), CoordRefSystems.ncoords(crs(grid)))
end

# ----------
# FALLBACKS
# ----------

Base.size(g::Grid) = size(topology(g))

paramdim(g::Grid) = length(size(g))

vertex(g::Grid, ind::Int) = vertex(g, CartesianIndices(vsize(g))[ind])

vertex(g::Grid, ijk::Dims) = vertex(g, LinearIndices(vsize(g))[ijk...])

Base.minimum(g::Grid) = vertex(g, ntuple(i -> 1, paramdim(g)))
Base.maximum(g::Grid) = vertex(g, vsize(g))
Base.extrema(g::Grid) = minimum(g), maximum(g)

function element(g::Grid, ind::Int)
  elem = element(topology(g), ind)
  type = pltype(elem)
  einds = indices(elem)
  cinds = CartesianIndices(vsize(g))
  verts = ntuple(i -> vertex(g, cinds[einds[i]]), nvertices(type))
  type(verts)
end

Base.eltype(g::Grid) = typeof(first(g))

Base.getindex(g::Grid, ind::Int) = element(g, ind)

Base.getindex(g::Grid, inds::AbstractVector) = [element(g, ind) for ind in inds]

Base.getindex(g::Grid, ijk::Int...) = element(g, LinearIndices(size(g))[ijk...])

@propagate_inbounds function Base.getindex(g::Grid, ijk...)
  dims = size(g)
  ranges = ntuple(i -> _asrange(dims[i], ijk[i]), paramdim(g))
  getindex(g, CartesianIndices(ranges))
end

function Base.getindex(g::Grid, I::CartesianIndices)
  @boundscheck _checkbounds(g, I)
  dims = size(I)
  odims = size(g)
  cinds = first(I):CartesianIndex(Tuple(last(I)) .+ 1)
  inds = vec(LinearIndices(odims .+ 1)[cinds])
  points = [vertex(g, ind) for ind in inds]
  periodic = isperiodic(topology(g)) .&& dims .== odims
  SimpleMesh(points, GridTopology(dims, periodic))
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

include("meshes/submeshes.jl")
include("meshes/regulargrid.jl")
include("meshes/cartesiangrid.jl")
include("meshes/rectilineargrid.jl")
include("meshes/structuredgrid.jl")
include("meshes/simplemesh.jl")
include("meshes/transfmeshes.jl")

# aliases for dispatch purposes
const OrthoRegularGrid{M<:𝔼,C<:Union{Cartesian,Projected}} = RegularGrid{M,C}
const OrthoRectilinearGrid{M<:𝔼,C<:Union{Cartesian,Projected}} = RectilinearGrid{M,C}
const OrthoStructuredGrid{M<:𝔼,C<:Union{Cartesian,Projected}} = StructuredGrid{M,C}
