# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Repair(K)

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
- K = 10: outer rings of polygon are expanded
- K = 11: rings of polygon are coherently oriented
- K = 12: degenerate rings of polygon are removed

## Examples

```
# remove duplicates and degenerates
mesh |> Repair(0) |> Repair(3)
```
"""
struct Repair{K} <: GeometricTransform end

Repair(K) = Repair{K}()

# --------------
# OPERATION (0)
# --------------

apply(::Repair{0}, geom::Polytope) = unique(geom), nothing

function apply(::Repair{0}, g::Chain)
  verts = vertices(g)

  # tuple of already seen edges (a, b) where a < b
  seen = Set{Tuple{eltype(verts),eltype(verts)}}()
  # already costructed chain parts
  parts = Vector{eltype(verts)}[]
  # currently working chain part
  current = [first(verts)]

  for i in 1:(length(verts) - 1)
    a = verts[i]
    b = verts[i + 1]

    # create a key for the edge that is independent of the order of vertices
    key = isless(a, b) ? (a, b) : (b, a)

    # if the edge has already been seen, we close the current part and start a new one
    if key in seen
      length(current) >= 2 && push!(parts, current)
      current = [b]
      # if the edge is new, we add it to the seen set and continue building the current part
    else
      push!(seen, key)
      push!(current, b)
    end
  end

  # if the last part has at least two vertices, we add it to the parts
  length(current) >= 2 && push!(parts, current)

  # construct geometries from the parts
  geoms = map(parts) do verts
    if length(verts) == 2
      Segment(verts...)
    elseif first(verts) == last(verts)
      Ring(verts[1:(end - 1)])
    else
      Rope(verts)
    end
  end

  maybemulti(geoms), nothing
end

function apply(::Repair{0}, g::PolyArea)
  repaired = map(rings(g)) do ring
    first(apply(Repair(0), ring))
  end

  outer = first(repaired)

  # Polygon collapsed to a lower-dimensional geometry
  outer isa Ring || return outer, nothing

  validholes = filter(r -> r isa Ring, repaired[2:end])

  PolyArea(outer, validholes), nothing
end

function apply(::Repair{0}, g::Ngon)
  repaired = first(apply(Repair(0), Ring(vertices(g))))

  repaired isa Ring || return repaired

  Ngon(vertices(repaired)...)
end

function _flatten!(out, geom)
    if geom isa Multi
        for g in parent(geom)
            _flatten!(out, g)
        end
    else
        push!(out, geom)
    end
end

function apply(::Repair{0}, g::Multi)
  repaired = map(parent(g)) do geom
    first(apply(Repair(0), geom))
  end

  flattened = Vector{eltype(repaired)}()
  for geom in repaired
    _flatten!(flattened, geom)
  end
  
  maybemulti(unique(flattened)), nothing
end

function apply(::Repair{0}, mesh::Mesh)
  # retrieve vertices and connectivities
  verts = vertices(mesh)
  elems = elements(topology(mesh))

  # remove duplicate vertices
  uverts = unique(verts)
  uindex = indexin(verts, uverts)

  # update indices of connectivities
  uelems = map(elems) do elem
    newinds = map(i -> uindex[i], indices(elem))
    connect(newinds, pltype(elem))
  end

  # remove duplicate connectivities
  uconnec = unique(elem -> Set(indices(elem)), uelems)

  # repaired mesh without duplicates
  rmesh = SimpleMesh(uverts, uconnec)

  rmesh, nothing
end

# --------------
# OPERATION (1)
# --------------

function apply(::Repair{1}, mesh::Mesh)
  # Ref to avoid boxing from map closure
  # see https://github.com/JuliaLang/julia/issues/15276
  count = Ref(0)

  # ordered list of indices
  seen = Int[]
  inds = Dict{Int,Int}()

  # reserve worst case memory need: all vertices are actually used
  sizehint!(seen, nvertices(mesh))
  sizehint!(inds, nvertices(mesh))

  topo = topology(mesh)
  elems = map(elements(topo)) do e
    elem = indices(e)
    for v in elem
      if !haskey(inds, v)
        push!(seen, v)
        count[] += 1
        inds[v] = count[]
      end
    end
    ntuple(i -> inds[elem[i]], length(elem))
  end

  points = [vertex(mesh, ind) for ind in seen]

  connec = map(connect, elems)

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

repair8(v::AbstractVector) = repair8(CircularVector(v))

function repair8(v::CircularVector{<:Point})
  n = length(v)
  keep = Int[]
  for i in 1:n
    t = Triangle(flat(v[i - 1]), flat(v[i]), flat(v[i + 1]))
    ispositive(area(t)) && push!(keep, i)
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

function _stretch10(g::Geometry)
  T = numtype(lentype(g))
  Stretch(ntuple(i -> one(T) + 10atol(T), embeddim(g)))
end

# ---------------
# OPERATION (11)
# ---------------

function apply(::Repair{11}, poly::PolyArea)
  r = rings(poly)

  # fix orientation
  ofix(r, o) = orientation(r) == o ? r : reverse(r)
  outer = ofix(first(r), CCW)
  inners = ofix.(r[2:end], CW)

  PolyArea([outer; inners]), nothing
end

# ---------------
# OPERATION (12)
# ---------------

function apply(::Repair{12}, poly::PolyArea)
  r = rings(poly)

  # fix degeneracy
  oring = first(r)
  outer = if nvertices(oring) == 2
    A, B = vertices(oring)
    P = centroid(Segment(A, B))
    Ring(A, P, B)
  else
    oring
  end

  # remove degenerated rings
  inners = filter(r -> nvertices(r) > 2, r[2:end])

  PolyArea([outer; inners]), nothing
end

# ----------
# FALLBACKS
# ----------

apply(::Repair, geom::Geometry) = geom, nothing

apply(t::Repair, multi::Multi) = Multi([t(g) for g in parent(multi)]), nothing

apply(t::Repair, dom::Domain) = GeometrySet([t(g) for g in dom]), nothing

# -----------
# IO METHODS
# -----------

Base.show(io::IO, ::Repair{K}) where {K} = print(io, "Repair(K: $K)")

function Base.show(io::IO, ::MIME"text/plain", t::Repair{K}) where {K}
  summary(io, t)
  println(io)
  print(io, "└─ K: $K")
end
