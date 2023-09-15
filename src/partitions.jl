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

Partition(object, subsets, metadata=Dict()) = Partition{typeof(object)}(object, subsets, metadata)

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

Base.iterate(partition::Partition, state=1) = state > length(partition) ? nothing : (partition[state], state + 1)

Base.length(partition::Partition) = length(partition.subsets)

Base.getindex(partition::Partition, ind::Int) = view(partition.object, partition.subsets[ind])

Base.getindex(partition::Partition, inds::AbstractVector{Int}) = [getindex(partition, ind) for ind in inds]

Base.eltype(partition::Partition) = typeof(first(partition))

function Base.summary(io::IO, partition::Partition)
  N = length(partition.subsets)
  print(io, "$N Partition")
end

Base.show(io::IO, partition::Partition) = summary(io, partition)

function Base.show(io::IO, ::MIME"text/plain", partition::Partition)
  meta = partition.metadata
  summary(io, partition)
  println(io)
  printelms(io, partition)
  if !isempty(meta)
    print(io, "\nmetadata: ")
    join(io, keys(meta), ", ")
  end
end
