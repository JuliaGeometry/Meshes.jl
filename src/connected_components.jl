# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    connected_components(mesh)

Return the connected components of `mesh` as a list of meshes.
"""
function connected_components(mesh::Mesh)
  graph = _meshgraph(mesh)
  components = connected_components(graph)
  [_meshcomponent(mesh, c) for c in components]
end

# construct the graph associated to a mesh
function _meshgraph(mesh)
  # get edges
  topo = convert(HalfEdgeTopology, topology(mesh))
  ∂₁₀ = Boundary{1,0}(topo)

  # make the graph
  graph = SimpleGraph(nvertices(mesh))
  for eind in 1:nfacets(topo)
    i, j = ∂₁₀(eind)
    add_edge!(graph, i, j)
  end

  graph
end

# make the mesh corresponding to a connected component
# (a component is a vector of vertex indices)
function _meshcomponent(mesh, component)
  # dictionary to map old indices to new indices
  newinds = Dict{Int,Int}()
  for (j, index) in enumerate(component)
    newinds[index] = j
  end

  # get faces
  topo = convert(HalfEdgeTopology, topology(mesh))
  ∂₂₀ = Boundary{2,0}(topo)

  # make the connectivities of the mesh component
  connec = Vector{Connectivity}(undef, 0)
  for faceind in 1:nelements(topo)
    faceinds = ∂₂₀(faceind) 
    if faceinds[1] in component
      facenewinds = [newinds[i] for i in faceinds]
      connectivity = connect(tuple(facenewinds...))
      push!(connec, connectivity)
    end
  end

  SimpleMesh(vertices(mesh)[component], connec)
end



