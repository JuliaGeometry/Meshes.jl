# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Partition(object, subsets, [metadata])

A partition of a spatial `object` into `subsets`.
Optionally, save `metadata` as a dictionary.
"""
struct Partition{O}
  object::O
  subsets::Vector{Vector{Int}}
  metadata::Dict
end

Partition(object, subsets, metadata=Dict()) =
  Partition{typeof(object)}(object, subsets, metadata)

"""
    indices(partition)

Return the subsets of indices in spatial object
that make up the `partition`.
"""
indices(partition::Partition) = partition.subsets

"""
    metadata(partition)

Return the metadata dictionary saved in the partition.
"""
metadata(partition::Partition) = partition.metadata

Base.iterate(partition::Partition, state=1) =
  state > length(partition) ? nothing : (partition[state], state + 1)

Base.length(partition::Partition) = length(partition.subsets)

Base.getindex(partition::Partition, ind::Int) =
  view(partition.object, partition.subsets[ind])

Base.getindex(partition::Partition, inds::AbstractVector{Int}) =
  [getindex(partition, ind) for ind in inds]

function Base.show(io::IO, partition::Partition)
  Dim = embeddim(partition.object)
  T   = coordtype(partition.object)
  N   = length(partition.subsets)
  print(io, "$N Partition{$Dim,$T}")
end

function Base.show(io::IO, ::MIME"text/plain", partition::Partition)
  subs = partition.subsets
  meta = partition.metadata
  println(io, partition)
  N = length(subs)
  I, J = N > 10 ? (5, N-4) : (N, N+1)
  lines = [["  └─$(partition[i])" for i in 1:I]
           (N > 10 ? ["  ⋮"] : [])
           ["  └─$(partition[i])" for i in J:N]]
  print(io, join(lines, "\n"))
  !isempty(meta) && print(io, "\n  metadata: ", join(keys(meta), ", "))
end
