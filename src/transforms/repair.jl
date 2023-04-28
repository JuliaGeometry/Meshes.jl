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
  nelem = nelements(topo)

  # indices of used vertices will be stored in this vector
  used = Vector{Int}(undef, 0)

  # initialize dictionary to map old indices to new indices
  newindices = Dict{Int,Int}()

  # to store the connectivities of the cleaned mesh
  connec = Vector{Connectivity}(undef, 0)

  # nused = length(used) will be incremented
  nused = 0

  # iterate over elements
  for e in 1:nelem
    face = ∂₂₀(e)
    for index in face
      if index ∉ used
        push!(used, index)
        nused = nused + 1
        newindices[index] = nused
      end
    end
    c = connect(tuple([newindices[i] for i in face]...))
    push!(connec, c)
  end

  # output mesh
  points = vertices(mesh)[used]
  SimpleMesh(points, connec)
end