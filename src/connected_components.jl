# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    connected_components(mesh)

Return the connected components of `mesh` as a list of meshes.
"""
function connected_components(mesh::Mesh)
  # get faces
  topo = convert(HalfEdgeTopology, topology(mesh))
  ∂₂₀ = Boundary{2,0}(topo)
  elemsind = collect(1:nelements(topo))

  # get vertices
  points = vertices(mesh)
  
  # get connected components
  graph = _meshgraph(mesh)
  components = connected_components(graph)
  meshes = Vector{SimpleMesh}(undef, length(components))

  # iterate over the connected components
  for (i, component) in enumerate(components)
    # dictionary to map old indices to new indices
    newinds = Dict{Int,Int}()
    for (j, index) in enumerate(component)
      newinds[index] = j
    end
    # vector to store the indices of the faces we will find, 
    # in order to delete them for the next iterations
    todelete = Vector{Int}(undef, 0)
    # make the connectivities of the mesh component
    connec = Vector{Connectivity}(undef, 0)
    for faceind in elemsind
      faceinds = ∂₂₀(faceind) 
      if faceinds[1] in component
        facenewinds = [newinds[i] for i in faceinds]
        connection = connect(tuple(facenewinds...))
        push!(connec, connection)
        push!(todelete, faceind)
      end
    end
    # make the mesh of this component
    meshes[i] = SimpleMesh(points[component], connec)
    # delete the indices of the faces we found
    setdiff!(elemsind, todelete)
  end

  meshes
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
