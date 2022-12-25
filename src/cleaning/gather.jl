# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    gather(mesh)

Merge the duplicated vertices of `mesh`.
"""
function gather(mesh)
  # get mesh vertices
  points = vertices(mesh)
  npoints = length(points)

  # vector indicating the duplicated vertices
  duplicated = fill(false, npoints)

  # dictionary to map the old indices to the new indices
  newindices = Dict{Int,Int}()

  # initialize new indices
  newindex = 1

  # iterate over vertices
  for index in 1:(npoints-1)
    if !duplicated[index]
      tail = points[(index+1):npoints]
      duplicates = index .+ findall(==(points[index]), tail)
      duplicated[duplicates] .= true
      push!(duplicates, index)
      for i in duplicates
        newindices[i] = newindex
      end
      newindex = newindex + 1
    end
  end

  # if the last vertex is not a duplicate, add it to the dictionary
  if !haskey(newindices, npoints)
    newindices[npoints] = newindex
  end

  # get the faces
  topo = convert(HalfEdgeTopology, topology(mesh))
  ∂₂₀ = Boundary{2,0}(topo)
  nfaces = nelements(topo)

  # vector to store the connectivities
  connectivities = Vector{Connectivity}(undef, 0)

  # iterate over the faces
  for f in 1:nfaces
    face = ∂₂₀(f)
    toconnect = [newindices[i] for i in face]
    connectivity = connect(tuple(toconnect...))
    push!(connectivities, connectivity)
  end

  SimpleMesh(points[.!duplicated], connectivities)
end
