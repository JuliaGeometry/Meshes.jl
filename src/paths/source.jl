# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SourcePath(sources, batchsize=10^3)

Traverse an object from `sources` and outwards.
Optionally pass a `batchsize` for KD-tree evaluations.
"""
struct SourcePath <: Path
  sources::Vector{Int}
  batchsize::Int
end

SourcePath(sources) = SourcePath(sources, 10^3)

function traverse(object, path::SourcePath)
  sources   = path.sources
  batchsize = path.batchsize
  @assert allunique(sources) "non-unique sources"
  @assert all(1 .≤ sources .≤ nelements(object)) "sources must be valid locations"
  @assert length(sources) ≤ nelements(object) "more sources than points in object"

  # coordinate matrix for source points
  S = coordinates(object, sources)

  # fit search tree
  kdtree = KDTree(S)

  # other locations that are not sources
  others = setdiff(1:nelements(object), sources)

  # process points in batches
  batches = Iterators.partition(others, batchsize)

  # pre-allocate memory for coordinates
  X = Matrix{coordtype(object)}(undef, embeddim(object), batchsize)

  # compute distances to sources
  dists = []
  for batch in batches
    coordinates!(X, object, batch)
    _, ds = knn(kdtree, view(X,:,1:length(batch)), length(sources), true)
    append!(dists, ds)
  end

  [sources; view(others, sortperm(dists))]
end
