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

# --------------
# OPERATION (1)
# --------------

function apply(::Repair{1}, mesh)
  count = 0
  seen = Int[]
  inds = Dict{Int,Int}()
  topo = topology(mesh)
  elems = map(elements(topo)) do e
    elem = indices(e)
    for v in elem
      if v ∉ seen
        push!(seen, v)
        count += 1
        inds[v] = count
      end
    end
    ntuple(i -> inds[elem[i]], length(elem))
  end

  points = vertices(mesh)[seen]

  connec = connect.(elems)

  rmesh = SimpleMesh(points, connec)

  rmesh, nothing
end

# --------------
# OPERATION (7)
# --------------

function apply(::Repair{7}, mesh)
  # HalfEdgeTopology constructor already
  # performs orientation of faces
  rmesh = topoconvert(HalfEdgeTopology, mesh)

  rmesh, nothing
end
