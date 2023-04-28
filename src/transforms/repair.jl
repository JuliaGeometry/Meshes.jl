# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Repair{K}

Perform repairing operation with code `K`.

## Available operations

- K = 0: duplicated vertices and faces are removed
- K = 1: unused vertices are removed
- K = 2: non-manifold faces are removed
- K = 3: degenerate faces are removed
- K = 4: non-manifold vertices are removed
- K = 5: non-manifold vertices are split by threshold
- K = 6: close vertices are merged (given a radius)
- K = 7: faces are coherently oriented
"""
struct Repair{K} <: StatelessGeometricTransform end

# implement operation K = 1
function apply(transform::Repair{1}, mesh)
  # get the faces
  topo = convert(HalfEdgeTopology, topology(mesh))
  ∂₂₀ = Boundary{2,0}(topo)
  nfaces = nelements(topo)

  # indices of used vertices will be stored in this vector
  used = Vector{Int}(undef, 0)

  # initialize dictionary to map old indices to new indices
  newindices = Dict{Int,Int}()

  # to store the connectivities of the cleaned mesh
  connectivities = Vector{Connectivity}(undef, 0)

  # l = length(used) will be incremented
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
    connec = connect(tuple([newindices[i] for i in face]...))
    push!(connectivities, connec)
  end

  # output mesh
  SimpleMesh(vertices(mesh)[used], connectivities)
end