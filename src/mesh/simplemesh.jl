# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SimpleMesh(vertices, topology)

A simple mesh with `vertices` and `topology`.

    SimpleMesh(vertices, connectivities; relations=false)

Alternatively, construct a simple mesh with `vertices` and
`connectivities`. The option `relations` can be used to
build topological relations assuming that the `connectivities`
represent the elements of the mesh.

## Examples

```julia
julia> points = [(0,0),(1,0),(1,1)]
julia> connec = [connect((1,2,3))]
julia> mesh   = SimpleMesh(points, connec)
```

See also [`Topology`](@ref), [`GridTopology`](@ref),
[`HalfEdgeTopology`](@ref), [`SimpleTopology`](@ref).

### Notes

- The option `relations=true` changes the underlying topology
  of the mesh to a [`HalfEdgeTopology`](@ref) instead of a
  [`SimpleTopology`](@ref).
"""
struct SimpleMesh{Dim,T,V<:AbstractVector{Point{Dim,T}},TP<:Topology} <: Mesh{Dim,T}
  vertices::V
  topology::TP
end

SimpleMesh(coords::AbstractVector{<:NTuple}, topology::Topology) =
  SimpleMesh(Point.(coords), topology)

function SimpleMesh(vertices, connec::AbstractVector{<:Connectivity}; relations=false)
  topology = relations ? HalfEdgeTopology(connec) : SimpleTopology(connec)
  SimpleMesh(vertices, topology)
end

vertices(m::SimpleMesh) = m.vertices

nvertices(m::SimpleMesh) = length(m.vertices)

"""
    convert(SimpleMesh, mesh)

Convert any `mesh` to a simple mesh with explicit
list of vertices and [`SimpleTopology`](@ref).
"""
Base.convert(::Type{<:SimpleMesh}, m::Mesh) = SimpleMesh(vertices(m), topology(m))
