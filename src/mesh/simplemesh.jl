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


"""
    isinside(point, mesh)

Tells whether a `point` is inside a surface `mesh`. Must be a surface mesh with only triangles as elements. Points that coincide with vertices of the mesh or lies on the face of a triangle return true.
"""
function isinside(testpoint::Point{3,T}, mesh::Mesh{3,T}) where T
  if !(eltype(mesh) <: Triangle)
    error("This function only works for surface meshes with triangles as elements.")
  end
  ex = testpoint .- extrema(mesh)
  direction = ex[argmax(norm.(ex))]
  r = Ray(testpoint, direction*2)
  
  intersects = false
  for elem in mesh
    if intersection(x -> x.type == NoIntersection ? false : true, r, elem)
      intersects = !intersects
    end
  end
  intersects
end
