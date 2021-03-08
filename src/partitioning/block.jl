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
  neighbors::Bool
end

BlockPartition(sides::NTuple{Dim,T}, neighbors = false) where {Dim,T} = BlockPartition{Dim,T}(sides, neighbors)
BlockPartition(sides::Vararg{T,Dim}; neighbors::Bool = false) where {Dim,T} = BlockPartition(SVector(sides), neighbors)

function partition(object, method::BlockPartition)
  Dim = embeddim(object)
  psides = method.sides
  bbox = boundingbox(object)
  calculate_metadata = method.neighbors
  
  println("Neighbors has been set to ", calculate_metadata)

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
  neighbors = [Int[] for i in 1:prod(nblocks)]

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

  if calculate_metadata == false
    println("Metadata calculation has been disabled. To calculate Metadata, please enable it while calling BlockPartition with neighbors = true.") 
  end
  #Intitialize metadata to an empty Dict.
  #If metadata calculation is enabled, we will populate the Dict.
  metadata = Dict()

  # neighboring blocks metadata
  if calculate_metadata == true
    println("Metadata Calculation has been enabled.")
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
    # Save Calculated Metadata
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
