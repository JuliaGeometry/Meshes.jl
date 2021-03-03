# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    BlockPartition(sides)
    BlockPartition(side₁, side₂, ...)

A method for partitioning spatial objects into blocks of given `sides`.
"""
struct BlockPartition{Dim,T} <: PartitionMethod
  sides::SVector{Dim,T}
end

BlockPartition(sides::NTuple{Dim,T}) where {Dim,T} = BlockPartition{Dim,T}(sides)

BlockPartition(sides::Vararg{T,Dim}) where {Dim,T} = BlockPartition(sides)

function partition(object, method::BlockPartition, calculate_metadata = false)
  Dim = embeddim(object)
  T = coordtype(object)

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

  subsets   = [Vector{Int}() for i in 1:prod(nblocks)]
  neighbors = [Vector{Int}() for i in 1:prod(nblocks)]

  # Cartesian to linear indices
  linear = LinearIndices(Dims(nblocks))

  coords = MVector{Dim,T}(undef)
  for j in 1:nelements(object)
    coordinates!(coords, object, j)

    # find block coordinates
    c = @. floor(Int, (coords - start) / psides) + 1
    @inbounds for i in 1:Dim
      c[i] = clamp(c[i], 1, nblocks[i])
    end
    bcoords = CartesianIndex(Tuple(c))

    # block index
    i = linear[bcoords]

    append!(subsets[i], j)
  end

  #Intitialize metadata to an empty Dict.
  #If calculate_metadata is enabled, we will populate it.
  metadata = Dict()

  # neighboring blocks metadata if calculate_metadata is enabled
  if calculate_metadata
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

    # save metadata if calculate_metadata is enabled
    metadata = Dict(:neighbors => neighbors)

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
