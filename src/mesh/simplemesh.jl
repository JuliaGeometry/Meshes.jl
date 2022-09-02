# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SimpleMesh(points, connec)

A simple mesh with `points` and connectivities `connec`.
The i-th face of the mesh is lazily built based on
the connectivity list `connec[i]`.

    SimpleMesh(points, topology)

Alternatively, construct a simple mesh with `points` and
a topological data structure (e.g. `HalfEdgeTopology`).

See also [`Topology`](@ref).

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
  SimpleMesh(points, FullTopology(connec))

vertices(m::SimpleMesh) = m.points

nvertices(m::SimpleMesh) = length(m.points)

topology(m::SimpleMesh) = m.topology

"""
    convert(SimpleMesh, mesh)

Convert any `mesh` to a simple mesh with explicit
list of points and [`FullTopology`](@ref).
"""
Base.convert(::Type{<:SimpleMesh}, m::Mesh) =
  SimpleMesh(vertices(m), topology(m))


"""
    point ∈ mesh

Tells whether a `point` is inside a surface `mesh`. Must be a surface mesh with only triangles as elements. Points that coincide with vertices of the mesh or lies on the face of a triangle return true.
"""
function Base.in(testpoint::Point{3,T}, mesh::Mesh{3,T}) where T
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
