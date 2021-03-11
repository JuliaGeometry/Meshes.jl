# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BlockPartition(sides; neighbors=false)
    BlockPartition(side₁, side₂, ...; neighbors=false)

A method for partitioning spatial objects into blocks of given `sides`.
Optionally, compute the `neighbors` of a block as the metadata.
"""
struct BlockPartition{Dim,T} <: PartitionMethod
  sides::SVector{Dim,T}
  neighbors::Bool
end

BlockPartition(sides::NTuple{Dim,T}; neighbors=false) where {Dim,T} = BlockPartition{Dim,T}(sides, neighbors)
BlockPartition(sides::Vararg{T,Dim}; neighbors=false) where {Dim,T} = BlockPartition(SVector(sides), neighbors)

function partition(object, method::BlockPartition)
  Dim = embeddim(object)

  psides = method.sides
  bbox = boundingbox(object)

  @assert all(psides .≤ sides(bbox)) "invalid block sides"

  # bounding box properties
  lo, up = coordinates.(extrema(bbox))
  ce = coordinates(center(bbox))

  # find number of blocks to left and right
  nleft  = @. ceil(Int, (ce - lo) / psides)
  nright = @. ceil(Int, (up - ce) / psides)

  start   = @. ce - nleft * psides
  nblocks = @. nleft + nright

  subsets   = [Int[] for i in 1:prod(nblocks)]

  # Cartesian to linear indices
  linear = LinearIndices(Dims(nblocks))

  for j in 1:nelements(object)
    coords = coordinates(centroid(object, j))

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
    bstart  = CartesianIndex(ntuple(i -> 1, Dim))
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

  Partition(object, subsets, metadata)
end
