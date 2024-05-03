# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BlockPartition(sides; neighbors=false)

A method for partitioning spatial objects into blocks of given `sides`.
Optionally, compute the `neighbors` of a block as the metadata.

    BlockPartition(side₁, side₂, ..., sideₙ; neighbors=false)

Alternatively, specify the sides `side₁`, `side₂`, ..., `sideₙ`.
"""
struct BlockPartition{Dim,ℒ<:Len} <: PartitionMethod
  sides::NTuple{Dim,ℒ}
  neighbors::Bool
  BlockPartition(sides::NTuple{Dim,ℒ}, neighbors::Bool) where {Dim,ℒ<:Len} = new{Dim,float(ℒ)}(sides, neighbors)
end

BlockPartition(sides::NTuple{Dim,Len}; neighbors=false) where {Dim} = BlockPartition(promote(sides...), neighbors)

BlockPartition(sides::Tuple; neighbors=false) = BlockPartition(addunit.(sides, u"m"), neighbors)

BlockPartition(sides...; neighbors=false) = BlockPartition(sides; neighbors=neighbors)

function partitioninds(::AbstractRNG, domain::Domain, method::BlockPartition)
  psides = method.sides

  bbox = boundingbox(domain)
  bsides = sides(bbox)
  Dim = length(bsides)

  @assert all(psides .≤ bsides) "invalid block sides"

  # bounding box properties
  ce = centroid(bbox)
  lo, up = extrema(bbox)

  # find number of blocks to left and right
  nleft = ceil.(Int, (ce - lo) ./ psides)
  nright = ceil.(Int, (up - ce) ./ psides)
  nblocks = @. nleft + nright

  # top left corner of first block
  start = coordinates(ce) .- nleft .* psides

  subsets = [Int[] for i in 1:prod(nblocks)]

  # Cartesian to linear indices
  linear = LinearIndices(Dims(nblocks))

  for j in 1:nelements(domain)
    coords = coordinates(centroid(domain, j))

    # find block coordinates
    c = @. floor(Int, (coords - start) / psides) + 1
    c = @. clamp(c, 1, nblocks)
    bcoords = CartesianIndex(Tuple(c))

    # block index
    i = linear[bcoords]

    append!(subsets[i], j)
  end

  # intitialize metadata
  metadata = Dict()
  neighbors = [Int[] for i in 1:prod(nblocks)]

  # neighboring blocks metadata
  if method.neighbors == true
    bstart = CartesianIndex(ntuple(i -> 1, Dim))
    boffset = CartesianIndex(ntuple(i -> 1, Dim))
    bfinish = CartesianIndex(Dims(nblocks))
    for (i, bcoords) in enumerate(bstart:bfinish)
      for b in (bcoords - boffset):(bcoords + boffset)
        if all(Tuple(bstart) .≤ Tuple(b) .≤ Tuple(bfinish)) && b ≠ bcoords
          push!(neighbors[i], linear[b])
        end
      end
    end

    metadata[:neighbors] = neighbors
  end

  # filter out empty blocks
  empty = isempty.(subsets)
  subsets = subsets[.!empty]
  neighbors = neighbors[.!empty]
  for i in findall(empty)
    for n in neighbors
      setdiff!(n, i)
    end
  end

  subsets, metadata
end
