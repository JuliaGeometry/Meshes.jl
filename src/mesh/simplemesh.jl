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


We can construct standard parametric surfaces from simpler meshes


```julia
function ϕ(u,v)

  x = 5/4*(1 - v/2π )*cos(2v)*(1 + cos(u)) + cos(2*v)
  y = 5/4*(1 - v/2π )*sin(2v)*(1 + cos(u)) + sin(2*v)
  z = 10*v/2π +5/4*(1 - v/2π )*sin(u) + 15
  return x,y,z

end

# Parametric Space
us = LinRange(0,2π,10)[1:end-1] #this component is periodic
vs = LinRange(0,π,40)

periodic = (false, true)
vertex = [(Point∘ϕ)(u,v)...) for u in us for v in vs]
topo = GridTopology( (length(vs)-1 , length(us)-1) .+ periodic, periodic)
mesh = SimpleMesh(vertex, topo)
viz(mesh, showfacets = true)
```

We can also use it to construct meshes from more complicated meshes

```julia
f(p) = p + rand(3) * 0.25
mesh_noisy = SimpleMesh([(Point∘f)(p) for p in coordinates.(vertices(mesh))], topology(mesh) )
viz(mesh_noisy)
```


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
  topology = relations ? HalfEdgeTopology(connec) :
                         SimpleTopology(connec)
  SimpleMesh(vertices, topology)
end

vertices(m::SimpleMesh) = m.vertices

nvertices(m::SimpleMesh) = length(m.vertices)

"""
    convert(SimpleMesh, mesh)

Convert any `mesh` to a simple mesh with explicit
list of vertices and [`SimpleTopology`](@ref).
"""
Base.convert(::Type{<:SimpleMesh}, m::Mesh) =
  SimpleMesh(vertices(m), topology(m))
