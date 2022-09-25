# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SimpleMesh(points, topology)

A simple mesh with `points` and `topology`.

    SimpleMesh(points, connec)

Alternatively, construct a simple mesh with `points` and
connectivities `connec` in a [`SimpleTopology`](@ref).

## Examples

```julia
julia> points = Point2[(0,0),(1,0),(1,1)]
julia> connec = connect.([(1,2,3)])
julia> SimpleMesh(points, connec)
```

See also [`Topology`](@ref), [`GridTopology`](@ref),
[`HalfEdgeTopology`](@ref), [`SimpleTopology`](@ref).

### Notes

- Connectivities must be given with coherent orientation, i.e.
  all faces must be counter-clockwise (CCW) or clockwise (CW).
"""
struct SimpleMesh{Dim,T,V<:AbstractVector{Point{Dim,T}},TP<:Topology} <: Mesh{Dim,T}
  points::V
  topology::TP
end

SimpleMesh(points::AbstractVector{<:Point},
           connec::AbstractVector{<:Connectivity}) =
  SimpleMesh(points, SimpleTopology(connec))

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
