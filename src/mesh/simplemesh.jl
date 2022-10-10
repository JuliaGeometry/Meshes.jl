# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SimpleMesh(points, topology)

A simple mesh with `points` and `topology`.

    SimpleMesh(points, connectivities; relations=false)

Alternatively, construct a simple mesh with `points` and
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
  points::V
  topology::TP
end

SimpleMesh(coords::AbstractVector{<:NTuple}, topology::Topology) =
  SimpleMesh(Point.(coords), topology)

function SimpleMesh(points, connec::AbstractVector{<:Connectivity}; relations=false)
  topology = relations ? HalfEdgeTopology(connec) :
                         SimpleTopology(connec)
  SimpleMesh(points, topology)
end

vertices(m::SimpleMesh) = m.points

nvertices(m::SimpleMesh) = length(m.points)

topology(m::SimpleMesh) = m.topology

"""
    convert(SimpleMesh, mesh)

Convert any `mesh` to a simple mesh with explicit
list of points and [`SimpleTopology`](@ref).
"""
Base.convert(::Type{<:SimpleMesh}, m::Mesh) =
  SimpleMesh(vertices(m), topology(m))
