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

function Partition(object, subsets, metadata=Dict())
  return Partition{typeof(object)}(object, subsets, metadata)
end

==(p₁::Partition, p₂::Partition) = p₁.object == p₂.object && p₁.subsets == p₂.subsets

Base.parent(p::Partition) = p.object

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

function Base.iterate(partition::Partition, state=1)
  return state > length(partition) ? nothing : (partition[state], state + 1)
end

Base.length(partition::Partition) = length(partition.subsets)

function Base.getindex(partition::Partition, ind::Int)
  return view(partition.object, partition.subsets[ind])
end

function Base.getindex(partition::Partition, inds::AbstractVector{Int})
  return [getindex(partition, ind) for ind in inds]
end

Base.eltype(partition::Partition) = typeof(first(partition))

function Base.show(io::IO, partition::Partition)
  N = length(partition.subsets)
  return print(io, "$N Partition")
end

function Base.show(io::IO, ::MIME"text/plain", partition::Partition)
  subs = partition.subsets
  meta = partition.metadata
  println(io, partition)
  N = length(subs)
  I, J = N > 10 ? (5, N - 4) : (N, N + 1)
  lines = [
    ["  └─$(partition[i])" for i in 1:I]
    (N > 10 ? ["  ⋮"] : [])
    ["  └─$(partition[i])" for i in J:N]
  ]
  print(io, join(lines, "\n"))
  return !isempty(meta) && print(io, "\n  metadata: ", join(keys(meta), ", "))
end
