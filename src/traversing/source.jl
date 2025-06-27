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
  assertion(allunique(sources), "non-unique sources")
  assertion(all(1 .≤ sources .≤ nelements(domain)), "sources must be valid locations")
  assertion(length(sources) ≤ nelements(domain), "more sources than points in object")

  # fit search tree
  kdtree = KDTree([svec(centroid(domain, s)) for s in sources])

  # other locations that are not sources
  others = setdiff(1:nelements(domain), sources)

  # process points in batches
  batches = Iterators.partition(others, batchsize)

  # compute distances to sources
  dists = map(batches) do batch
    xs = [svec(centroid(domain, b)) for b in batch]
    _, ds = nn(kdtree, xs)
    ds
  end

  # sort other indices
  alldists = reduce(vcat, dists)
  perminds = sortperm(alldists)
  permuted = view(others, perminds)

  [sources; permuted]
end
