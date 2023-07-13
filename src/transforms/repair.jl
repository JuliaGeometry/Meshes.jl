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
- K = 8: zero-area ears are removed
"""
struct Repair{K} <: StatelessGeometricTransform end

# --------------
# OPERATION (0)
# --------------

apply(::Repair{0}, geom::Polytope) = unique(geom), nothing

apply(::Repair{0}, mesh::Mesh) = @error "not implemented"

# --------------
# OPERATION (1)
# --------------

function apply(::Repair{1}, mesh::Mesh)
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

# HalfEdgeTopology constructor
# performs orientation of faces
apply(::Repair{7}, mesh::Mesh) = topoconvert(HalfEdgeTopology, mesh), nothing

# --------------
# OPERATION (8)
# --------------

function apply(::Repair{8}, poly::PolyArea)
  v = poly |> rings .|> vertices .|> repair8
  p = if hasholes(poly)
    PolyArea(v[begin], v[(begin + 1):end])
  else
    PolyArea(v[begin])
  end
  p, nothing
end

function apply(::Repair{8}, poly::Ngon)
  v = poly |> vertices |> repair8
  N = length(v)
  Ngon(ntuple(i -> @inbounds(v[i]), N)), nothing
end

function apply(::Repair{8}, ring::Ring)
  v = ring |> vertices |> repair8
  Ring(v), nothing
end

repair8(v) = repair8(collect(v))

repair8(v::AbstractVector) = repair8(CircularVector(v))

function repair8(v::CircularVector{Point{Dim,T}}) where {Dim,T}
  n = length(v)
  keep = Int[]
  for i in 1:n
    t = Triangle(v[i - 1], v[i], v[i + 1])
    area(t) > atol(T)^2 && push!(keep, i)
  end
  v[keep]
end
