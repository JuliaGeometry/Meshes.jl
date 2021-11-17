# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CartesianGrid(dims, reference, spacing)

A Cartesian grid with dimensions `dims`, lower left corner at `reference`
and cell spacing `spacing`. The three arguments must have the same length.

    CartesianGrid(start, finish, dims=dims)

Alternatively, construct a Cartesian grid from a `start` point (lower left)
to a `finish` point (upper right).

    CartesianGrid(start, finish, spacing)

Alternatively, construct a Cartesian grid from a `start` point to a `finish`
point using a given `spacing`.

    CartesianGrid(dims)
    CartesianGrid(dim1, dim2, ...)

Finally, a Cartesian grid can be constructed by only passing the dimensions
`dims` as a tuple, or by passing each dimension `dim1`, `dim2`, ... separately.
In this case, the reference and spacing default to (0,0,...) and (1,1,...).

## Examples

Create a 3D grid with 100x100x50 locations:

```julia
julia> CartesianGrid(100,100,50)
```

Create a 2D grid with 100x100 locations and reference at (10.,20.) units:

```julia
julia> CartesianGrid((100,100),(10.,20.),(1.,1.))
```

Create a 1D grid from -1 to 1 with 100 locations:

```julia
julia> CartesianGrid((-1.,),(1.,), dims=(100,))
```
"""
struct CartesianGrid{Dim,T} <: Mesh{Dim,T}
  dims::Dims{Dim}
  reference::Point{Dim,T}
  spacing::SVector{Dim,T}
  offset::SVector{Dim,Int}

  function CartesianGrid{Dim,T}(dims, reference, spacing, offset) where {Dim,T}
    @assert all(dims .> 0) "dimensions must be positive"
    @assert all(spacing .> 0) "spacing must be positive"
    new(dims, reference, spacing, offset)
  end
end

CartesianGrid(dims::Dims{Dim}, reference::Point{Dim,T},
              spacing::SVector{Dim,T},
              offset::SVector{Dim,Int}=SVector(ntuple(i->1, Dim))) where {Dim,T} =
  CartesianGrid{Dim,T}(dims, reference, spacing, offset)

CartesianGrid(dims::Dims{Dim}, reference::NTuple{Dim,T},
              spacing::NTuple{Dim,T},
              offset::NTuple{Dim,Int}=ntuple(i->1, Dim)) where {Dim,T} =
  CartesianGrid{Dim,T}(dims, Point(reference), SVector(spacing), SVector(offset))

function CartesianGrid(start::Point{Dim,T}, finish::Point{Dim,T},
                       spacing::SVector{Dim,T}) where {Dim,T}
  dims = Tuple(ceil.(Int, (finish - start) ./ spacing))
  CartesianGrid{Dim,T}(dims, start, spacing, SVector(ntuple(i->1, Dim)))
end

CartesianGrid(start::NTuple{Dim,T}, finish::NTuple{Dim,T},
              spacing::NTuple{Dim,T}) where {Dim,T} =
  CartesianGrid(Point(start), Point(finish), SVector(spacing))

CartesianGrid(start::Point{Dim,T}, finish::Point{Dim,T};
              dims::Dims{Dim}=ntuple(i->100, Dim)) where {Dim,T} =
  CartesianGrid{Dim,T}(dims, start, (finish - start) ./ dims, SVector(ntuple(i->1, Dim)))

CartesianGrid(start::NTuple{Dim,T}, finish::NTuple{Dim,T};
              dims::Dims{Dim}=ntuple(i->100, Dim)) where {Dim,T} =
  CartesianGrid(Point(start), Point(finish); dims=dims)

CartesianGrid{T}(dims::Dims{Dim}) where {Dim,T} =
  CartesianGrid{Dim,T}(dims, ntuple(i->zero(T), Dim), ntuple(i->one(T), Dim), ntuple(i->1, Dim))

CartesianGrid{T}(dims::Vararg{Int,Dim}) where {Dim,T} = CartesianGrid{T}(dims)

CartesianGrid(dims::Dims{Dim}) where {Dim} = CartesianGrid{Float64}(dims)

CartesianGrid(dims::Vararg{Int,Dim}) where {Dim} = CartesianGrid{Float64}(dims)

==(g1::CartesianGrid, g2::CartesianGrid) =
  g1.dims    == g2.dims    &&
  g1.reference  == g2.reference  &&
  g1.spacing == g2.spacing &&
  g1.offset == g2.offset

Base.size(g::CartesianGrid) = g.dims
Base.minimum(g::CartesianGrid) = g.reference - (g.offset .- 1) .* g.spacing
Base.maximum(g::CartesianGrid) = g.reference + (g.dims .+ 1 .- g.offset) .* g.spacing
Base.extrema(g::CartesianGrid) = minimum(g), maximum(g)
spacing(g::CartesianGrid) = g.spacing

function vertices(g::CartesianGrid)
  inds_ranges = map((dims,offset) -> (0:dims) .+ 1 .- offset, g.dims, g.offset.data)
  inds = CartesianIndices(inds_ranges)
  vec([g.reference + ind.I .* g.spacing for ind in inds])
end

elements(g::CartesianGrid) = (g[i] for i in 1:nelements(g))

topology(g::CartesianGrid) = GridTopology(size(g))

# -----------------
# DOMAIN INTERFACE
# -----------------

function element(g::CartesianGrid{Dim}, ind::Int) where {Dim}
  inds_ranges = map((dims,offset) -> (1:dims) .+ 1 .- offset, g.dims, g.offset.data)
  I = CartesianIndices(inds_ranges)[ind]
  o = coordinates(g.reference)
  s = g.spacing
  i = I.I

  if Dim == 1 # segment
    p1 = (o[1] + (i[1] - 1) * s[1],)
    p2 = (o[1] + (    i[1]) * s[1],)
    Segment(p1, p2)
  elseif Dim == 2 # quadrangle
    p1 = (o[1] + (i[1] - 1) * s[1],
          o[2] + (i[2] - 1) * s[2])
    p2 = (o[1] + (i[1]    ) * s[1],
          o[2] + (i[2] - 1) * s[2])
    p3 = (o[1] + (i[1]    ) * s[1],
          o[2] + (i[2]    ) * s[2])
    p4 = (o[1] + (i[1] - 1) * s[1],
          o[2] + (i[2]    ) * s[2])
    Quadrangle(p1, p2, p3, p4)
  elseif Dim == 3 # hexahedron
    p1 = (o[1] + (i[1] - 1) * s[1],
          o[2] + (i[2] - 1) * s[2],
          o[3] + (i[3] - 1) * s[3])
    p2 = (o[1] + (i[1]    ) * s[1],
          o[2] + (i[2] - 1) * s[2],
          o[3] + (i[3] - 1) * s[3])
    p3 = (o[1] + (i[1]    ) * s[1],
          o[2] + (i[2]    ) * s[2],
          o[3] + (i[3] - 1) * s[3])
    p4 = (o[1] + (i[1] - 1) * s[1],
          o[2] + (i[2]    ) * s[2],
          o[3] + (i[3] - 1) * s[3])
    p5 = (o[1] + (i[1] - 1) * s[1],
          o[2] + (i[2] - 1) * s[2],
          o[3] + (i[3]    ) * s[3])
    p6 = (o[1] + (i[1]    ) * s[1],
          o[2] + (i[2] - 1) * s[2],
          o[3] + (i[3]    ) * s[3])
    p7 = (o[1] + (i[1]    ) * s[1],
          o[2] + (i[2]    ) * s[2],
          o[3] + (i[3]    ) * s[3])
    p8 = (o[1] + (i[1] - 1) * s[1],
          o[2] + (i[2]    ) * s[2],
          o[3] + (i[3]    ) * s[3])
    Hexahedron(p1, p2, p3, p4, p5, p6, p7, p8)
  else
    throw(ErrorException("not implemented"))
  end
end

nelements(g::CartesianGrid) = prod(g.dims)

function centroid(g::CartesianGrid{Dim}, ind::Int) where {Dim}
  inds_ranges = map((dims,offset) -> (1:dims) .+ 1 .- offset, g.dims, g.offset.data)
  intcoords = CartesianIndices(inds_ranges)[ind]
  newreference = coordinates(g.reference) .+ g.spacing ./ 2
  Point(ntuple(i -> newreference[i] + (intcoords[i] - 1)*g.spacing[i], Dim))
end

Base.eltype(g::CartesianGrid) = typeof(g[1])

# ----------------
# OPTIONAL TRAITS
# ----------------

isgrid(::Type{<:CartesianGrid}) = true

# ----------------------------
# ADDITIONAL INDEXING METHODS
# ----------------------------

"""
    grid[istart:iend,jstart:jend,...]

Return a subgrid of the Cartesian `grid` using integer ranges
`istart:iend`, `jstart:jend`, ...
"""
Base.getindex(g::CartesianGrid{Dim}, r::Vararg{UnitRange{Int},Dim}) where {Dim} =
  getindex(g, CartesianIndex(first.(r)):CartesianIndex(last.(r)))

function Base.getindex(g::CartesianGrid{Dim}, I::CartesianIndices{Dim}) where {Dim}
  offset = g.offset .- first(I).I .+ 1
  dims   = size(I)
  CartesianGrid(dims, g.reference, g.spacing, offset)
end

Base.view(g::CartesianGrid{Dim}, I::CartesianIndices{Dim}) where {Dim} = getindex(g, I)

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, g::CartesianGrid{Dim,T}) where {Dim,T}
  dims = join(g.dims, "×")
  print(io, "$dims CartesianGrid{$Dim,$T}")
end

function Base.show(io::IO, ::MIME"text/plain", g::CartesianGrid)
  println(io, g)
  println(io, "  minimum: ", minimum(g))
  println(io, "  maximum: ", maximum(g))
  print(  io, "  spacing: ", Tuple(spacing(g)))
end
