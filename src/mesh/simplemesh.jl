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
julia> points = [(0.0, 0.0),(1.0, 0.0), (1.0, 1.0)]
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
struct SimpleMesh{Dim,V<:AbstractVector{<:Point{Dim}},TP<:Topology} <: Mesh{Dim,TP}
  vertices::V
  topology::TP
end

SimpleMesh(coords::AbstractVector{<:Tuple}, topology::Topology) = SimpleMesh(Point.(coords), topology)

function SimpleMesh(vertices, connec::AbstractVector{<:Connectivity}; relations=false)
  topology = relations ? HalfEdgeTopology(connec) : SimpleTopology(connec)
  SimpleMesh(vertices, topology)
end

lentype(::Type{<:SimpleMesh{Dim,V}}) where {Dim,V} = lentype(eltype(V))

vertex(m::SimpleMesh, ind::Int) = m.vertices[ind]

vertices(m::SimpleMesh) = m.vertices

nvertices(m::SimpleMesh) = length(m.vertices)

function Base.getindex(m::SimpleMesh{Dim,V,GridTopology{Dim}}, I::CartesianIndices{Dim}) where {Dim,V}
  @boundscheck _checkbounds(m, I)
  dims = size(I)
  odims = size(m)
  cinds = first(I):CartesianIndex(Tuple(last(I)) .+ 1)
  inds = vec(LinearIndices(odims .+ 1)[cinds])
  periodic = isperiodic(topology(m)) .&& dims .== odims
  SimpleMesh(m.vertices[inds], GridTopology(dims, periodic))
end
