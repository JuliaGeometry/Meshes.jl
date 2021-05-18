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
  println(io, _lines(verts, "    "))
  println(io, "  $nelms elements")
  print(  io, _lines(elems, "    "))
end

function _lines(itr, tab="  ")
  vec = collect(itr)
  N = length(vec)
  I, J = N > 10 ? (5, N-4) : (N, N+1)
  lines = [["$(tab)└─$(vec[i])" for i in 1:I]
           (N > 10 ? ["$(tab)⋮"] : [])
           ["$(tab)└─$(vec[i])" for i in J:N]]
  join(lines, "\n")
end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("mesh/cartesiangrid.jl")
include("mesh/simplemesh.jl")
