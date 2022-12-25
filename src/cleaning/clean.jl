# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    clean(mesh)

Remove the unused vertices of `mesh`.
"""
function clean(mesh)
  # get the faces of the mesh
  topo = convert(HalfEdgeTopology, topology(mesh))
  ∂₂₀ = Boundary{2,0}(topo)
  nfaces = nelements(topo)

  # vector to store the indices of the used vertices
  used = Vector{Int}(undef, 0)

  # dictionary to map the old indices to the new indices
  newindices = Dict{Int,Int}()

  # vector to store the connectivities of the cleaned mesh
  connectivities = Vector{Connectivity}(undef, 0)

  # l=length(used) will be incremented
  l = 0

  # iterate over faces
  for f in 1:nfaces
    face = ∂₂₀(f)
    for index in face
      if !in(index, used)
        push!(used, index)
        l = l + 1
        newindices[index] = l
      end
    end
    connectivity = connect(tuple([newindices[i] for i in face]...))
    push!(connectivities, connectivity)
  end

  SimpleMesh(vertices(mesh)[used], connectivities)
end
