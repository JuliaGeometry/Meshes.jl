# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Mesh{Dim,T}

A mesh embedded in a `Dim`-dimensional space with coordinates of type `T`.
"""
abstract type Mesh{Dim,T} <: Domain{Dim,T} end

"""
    topology(mesh)

Return the topological structure of the `mesh`.
"""
function topology end

"""
    vertices(mesh)

Return the vertices of the `mesh`.
"""
function vertices end

"""
    nvertices(mesh)

Return the number of vertices of the `mesh`.
"""
nvertices(m::Mesh) = nvertices(topology(m))

"""
    faces(mesh, rank)

Return an iterator with `rank`-faces of the `mesh`.

## Example

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
    elements(mesh)

Return the top-faces (a.k.a. elements) of the `mesh`.

## Example

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

Return the number of elements in the `mesh`.
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

Return the number of facets in the `mesh`.
"""
nfacets(m::Mesh) = nfacets(topology(m))

"""
    topoconvert(T, mesh)

Convert underlying topology of the `mesh` to topology of type `T`.

## Example

Convert underlying topology to [`HalfEdgeTopology`](@ref) for
efficient topological relations.

```julia
newmesh = topoconvert(HalfEdgeTopology, mesh)
```
"""
topoconvert(TP::Type{<:Topology}, m::Mesh) =
  SimpleMesh(vertices(m), convert(TP, topology(m)))

"""
    merge(mesh₁, mesh₂)

Merge `mesh₁` with `mesh₂`, i.e. concatenate the vertices
and adjust the connectivities accordingly.
"""
function Base.merge(m₁::Mesh, m₂::Mesh)
  t₁ = topology(m₁)
  t₂ = topology(m₂)
  v₁ = collect(vertices(m₁))
  v₂ = collect(vertices(m₂))
  e₁ = collect(elements(t₁))
  e₂ = collect(elements(t₂))

  # concatenate vertices
  points = [v₁; v₂]

  # concatenate connectivities
  connec = e₁
  offset = length(v₁)
  for e in e₂
    c  = indices(e)
    c′ = ntuple(i -> c[i] + offset, length(c))
    push!(connec, connect(c′))
  end

  SimpleMesh(points, connec)
end

==(m₁::Mesh, m₂::Mesh) =
  vertices(m₁) == vertices(m₂) &&
  topology(m₁) == topology(m₂)

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, ::MIME"text/plain", m::Mesh{Dim,T}) where {Dim,T}
  t = topology(m)
  verts = vertices(m)
  elems = elements(t)
  nvert = nvertices(m)
  nelms = nelements(m)
  println(io, m)
  println(io, "  $nvert vertices")
  println(io, io_lines(verts, "    "))
  println(io, "  $nelms elements")
  print(  io, io_lines(elems, "    "))
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("mesh/cartesiangrid.jl")
include("mesh/simplemesh.jl")
