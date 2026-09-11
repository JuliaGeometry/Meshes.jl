# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    constructor(G)

Given a (parametric) type `G{T₁,T₂,...}`, return the type `G`.
"""
constructor(G::Type) = getfield(Meshes, nameof(G))

"""
    fitdims(dims, D)

Fit tuple `dims` to a given length `D` by repeating the last dimension.
"""
function fitdims(dims::Dims{N}, D) where {N}
  ntuple(i -> i ≤ N ? dims[i] : last(dims), D)
end

"""
    collectat(iter, inds)

Collect iterable `iter` at indices `inds` efficiently.
"""
function collectat(iter, inds)
  if isempty(inds)
    eltype(iter)[]
  else
    m = maximum(inds)
    e = Iterators.enumerate(iter)
    w = Iterators.takewhile(x -> (first(x) ≤ m), e)
    f = Iterators.filter(x -> (first(x) ∈ inds), w)
    map(last, f)
  end
end

collectat(vec::AbstractVector, inds) = vec[inds]

"""
    XYZ(xyz)

Generate the coordinate arrays `XYZ` from the coordinate vectors `xyz`.
"""
@generated function XYZ(xyz::NTuple{Dim,AbstractVector}) where {Dim}
  exprs = ntuple(Dim) do d
    quote
      a = xyz[$d]
      A = Array{eltype(a),Dim}(undef, length.(xyz))
      @nloops $Dim i A begin
        @nref($Dim, A, i) = a[$(Symbol(:i_, d))]
      end
      A
    end
  end
  Expr(:tuple, exprs...)
end

"""
    isthreaded(cond=true)

Return true if `cond`ition is true in the presence of multiple threads.
"""
isthreaded(cond=true) = cond && Threads.nthreads() > 1

"""
    maybemulti(geoms)

Return `only(geoms)` or `Multi(geoms)` depending on the length of the `geoms` vector.
"""
maybemulti(geoms::AbstractVector{<:Geometry}) = length(geoms) == 1 ? only(geoms) : Multi(geoms)

"""
    glue(segs)

Glue unique segments into `Segment`s, `Rope`s and `Ring`s. Segments are sorted.
"""
function glue(segs::AbstractVector{<:Segment})
  length(segs) == 1 && return [only(segs)]

  # sort segments independently of input order
  segs = sort(segs; by=seg -> begin
    a, b = vertices(seg)
    ifelse(a < b, (a, b), (b, a))
  end)

  # determine the manifold and crs of the segments
  M = manifold(eltype(segs))
  C = crs(eltype(segs))

  # build adjacency dictionary
  adj = Dict{Point{M,C},Vector{Int}}()
  for (i, seg) in enumerate(segs)
    a, b = vertices(seg)
    push!(get!(adj, a, Int[]), i)
    push!(get!(adj, b, Int[]), i)
  end

  visited = falses(length(segs))
  chains = Chain{M,C}[]

  # trace a maximal path from a starting vertex through a segment
  function trace(start, segind)
    verts = [start]
    current = start
    currentind = segind

    while true
      visited[currentind] = true

      a, b = vertices(segs[currentind])
      next = ifelse(current ≈ a, b, a)
      push!(verts, next)

      # stop at terminal or branching vertices
      length(adj[next]) == 2 || break

      i, j = adj[next]
      nextind = ifelse(i == currentind, j, i)

      visited[nextind] && break

      current = next
      currentind = nextind
    end

    verts
  end

  # first trace maximal non-cyclic paths
  starts = sort([v for (v, inds) in adj if length(inds) != 2])

  for start in starts
    for segind in adj[start]
      visited[segind] && continue

      verts = trace(start, segind)

      geom = length(verts) == 2 ? Segment(verts...) : Rope(verts)
      push!(chains, geom)
    end
  end

  # remaining unvisited segments belong to cycles
  for segind in eachindex(segs)
    visited[segind] && continue

    a, b = vertices(segs[segind])
    start = ifelse(a < b, a, b)

    verts = trace(start, segind)

    push!(chains, Ring(verts[1:(end - 1)]))
  end

  chains
end

function glue(g::Multi)
  geoms = flatten(g)

  # separate between points and chains
  points = unique(filter(geom -> geom isa Point, geoms))
  chains = filter(geom -> geom isa Chain, geoms)

  # if there are no chains, we can return the glued geometry immediately
  isempty(chains) && return points

  # collect all segments from the chains
  segs = collect(Iterators.flatten(segments(g) for g in chains))
  glued = glue(segs)

  # if there are no points, we can return the glued geometry immediately
  isempty(points) && return glued

  # remove points already represented by the glued 1D geometry
  points = filter(p -> !any(chain -> p ∈ chain, glued), points)

  [glued..., points...]
end

"""
    flatten(g)

Return a vector of the geometries contained in `g`. If `g` is a `Multi`, it recursively flattens its parents.
The output is a vector containing all the geometries within `g`, with any nested `Multi` geometries recursively flattened.

See [`maybemulti`](@ref) for turning this vector into a single geometry. 
"""
flatten(g::Multi) = mapreduce(flatten, vcat, parent(g))

flatten(g::Geometry) = [g]
