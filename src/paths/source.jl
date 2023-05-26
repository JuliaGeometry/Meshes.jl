# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SourcePath(sources, batchsize=10^3)

Traverse a domain from `sources` and outwards.
Optionally pass a `batchsize` for KD-tree evaluations.
"""
struct SourcePath <: Path
  sources::Vector{Int}
  batchsize::Int
end

SourcePath(sources) = SourcePath(sources, 10^3)

function traverse(domain, path::SourcePath)
  sources = path.sources
  batchsize = path.batchsize
  @assert allunique(sources) "non-unique sources"
  @assert all(1 .≤ sources .≤ nelements(domain)) "sources must be valid locations"
  @assert length(sources) ≤ nelements(domain) "more sources than points in object"

  # fit search tree
  kdtree = KDTree(coordinates.([centroid(domain, s) for s in sources]))

  # other locations that are not sources
  others = setdiff(1:nelements(domain), sources)

  # process points in batches
  batches = Iterators.partition(others, batchsize)

  # compute distances to sources
  dists = []
  for batch in batches
    coords = coordinates.([centroid(domain, b) for b in batch])
    _, ds = knn(kdtree, coords, length(sources), true)
    append!(dists, ds)
  end

  [sources; view(others, sortperm(dists))]
end
