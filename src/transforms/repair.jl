# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Repair{K}

Perform repairing operation with code `K`.

## Available operations

- K =  0: duplicated vertices and faces are removed
- K =  1: unused vertices are removed
- K =  2: non-manifold faces are removed
- K =  3: degenerate faces are removed
- K =  4: non-manifold vertices are removed
- K =  5: non-manifold vertices are split by threshold
- K =  6: close vertices are merged (given a radius)
- K =  7: faces are coherently oriented
- K =  8: zero-area ears are removed
- K =  9: rings of polygon are sorted
- K = 10: outer rings are expanded

## Examples

```
# remove duplicates and degenerates
mesh |> Repair{0}() |> Repair{3}()
```
"""
struct Repair{K} <: GeometricTransform end

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

# HalfEdgeTopology constructor performs orientation of faces
apply(::Repair{7}, mesh::Mesh) = topoconvert(HalfEdgeTopology, mesh), nothing

# --------------
# OPERATION (8)
# --------------

function apply(::Repair{8}, poly::PolyArea)
  v = poly |> rings .|> vertices .|> repair8
  PolyArea(v), nothing
end

function apply(::Repair{8}, poly::Ngon)
  v = poly |> vertices |> repair8
  Ngon(ntuple(i -> @inbounds(v[i]), length(v))), nothing
end

function apply(::Repair{8}, ring::Ring)
  v = ring |> vertices |> repair8
  Ring(v), nothing
end

repair8(v) = repair8(collect(v))

repair8(v::AbstractVector) = repair8(CircularVector(v))

function repair8(v::CircularVector{<:Point})
  n = length(v)
  keep = Int[]
  for i in 1:n
    t = Triangle(v[i - 1], v[i], v[i + 1])
    a = area(t)
    a > atol(a) && push!(keep, i)
  end
  isempty(keep) ? v[begin] : v[keep]
end

# --------------
# OPERATION (9)
# --------------

function apply(::Repair{9}, poly::PolyArea)
  newrings, indices = poly |> rings |> repair9
  PolyArea(newrings), indices
end

apply(::Repair{9}, poly::Ngon) = poly, []

function repair9(r::AbstractVector{<:Ring})
  # sort vertices lexicographically
  verts = vertices.(r)
  coord = to.(reduce(vcat, verts))
  vperm = sortperm(sortperm(coord))

  # each ring has its own set of indices
  offset = 0
  indices = Vector{Int}[]
  for vert in verts
    nvert = length(vert)
    range = (offset + 1):(offset + nvert)
    push!(indices, vperm[range])
    offset += nvert
  end

  # sort rings based on leftmost vertex
  leftmost = argmin.(indices)
  minimums = getindex.(indices, leftmost)
  neworder = sortperm(minimums)
  newverts = verts[neworder]
  newinds = indices[neworder]

  Ring.(newverts), newinds
end

# ---------------
# OPERATION (10)
# ---------------

function apply(::Repair{10}, poly::PolyArea)
  t = _stretch10(poly)
  r = rings(poly)
  n, c = apply(t, first(r))
  PolyArea([n; r[2:end]]), (t, c)
end

function revert(::Repair{10}, poly::PolyArea, c)
  r = rings(poly)
  o = revert(c[1], first(r), c[2])
  PolyArea([o; r[2:end]])
end

apply(::Repair{10}, poly::Ngon) = poly, nothing

revert(::Repair{10}, poly::Ngon, cache) = poly

function _stretch10(g::Geometry)
  T = numtype(lentype(g))
  Stretch(ntuple(i -> one(T) + 10atol(T), embeddim(g)))
end
