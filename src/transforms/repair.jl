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
# OPERATION (0)
# --------------

function apply(::Repair{0}, mesh)
  points  = vertices(mesh)
  npoints = nvertices(mesh)

  count = 0
  dups = fill(false, npoints)
  inds = Dict{Int,Int}()

  for i in 1:(npoints-1)
    if !dups[i]
      # find duplicates
      js = i .+ findall(==(points[i]), points[i+1:npoints])
      dups[js] .= true
      count += 1
      for k in [js; i]
        inds[k] = count
      end
    end
  end
  
  # if last vertex is not a duplicate, add it to the dictionary
  if npoints ∉ keys(inds) 
    inds[npoints] = count
  end

  # get the elements (faces)
  topo = topology(mesh)
  ∂₂₀ = Boundary{2,0}(topo)
  nelems = nelements(topo)

  # vector to store the connectivities
  connec = Vector{Connectivity}(undef, 0)

  # iterate over the elements
  for e in 1:nelems
    elem = ∂₂₀(e)
    toconnect = [inds[i] for i in elem]
    c = connect(ntuple(i -> toconnect[i], length(toconnect)))
    push!(connec, c)
  end

  # unique points
  upoints = points[.!dups]

  # indices of non-duplicates faces
  uinds = unique(i -> Set{Int}([indices(connec[i])...]), 1:nelems)

  # unique connectivities
  uconnec = connec[uinds]

  rmesh = SimpleMesh(upoints, uconnec)

  rmesh, nothing
end

# --------------
# OPERATION (1)
# --------------

function apply(::Repair{1}, mesh)
  topo = convert(HalfEdgeTopology, topology(mesh))
  ∂₂₀  = Boundary{2,0}(topo)

  count = 0
  seen  = Int[]
  inds  = Dict{Int,Int}()
  elems = map(1:nelements(mesh)) do e
    elem = ∂₂₀(e)
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