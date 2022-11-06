# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CartesianGrid(dims, origin, spacing)

A Cartesian grid with dimensions `dims`, lower left corner at `origin`
and cell spacing `spacing`. The three arguments must have the same length.

    CartesianGrid(dims, origin, spacing, offset)

A Cartesian grid with dimensions `dims`, with lower left corner of element
`offset` at `origin` and cell spacing `spacing`.

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
In this case, the origin and spacing default to (0,0,...) and (1,1,...).

## Examples

Create a 3D grid with 100x100x50 hexahedrons:

```julia
julia> CartesianGrid(100,100,50)
```

Create a 2D grid with 100x100 quadrangles and origin at (10.,20.) units:

```julia
julia> CartesianGrid((100,100),(10.,20.),(1.,1.))
```

Create a 1D grid from -1 to 1 with 100 segments:

```julia
julia> CartesianGrid((-1.,),(1.,), dims=(100,))
```
"""
struct CartesianGrid{Dim,T} <: Mesh{Dim,T}
  origin::Point{Dim,T}
  spacing::NTuple{Dim,T}
  offset::Dims{Dim}
  topology::GridTopology{Dim}

  function CartesianGrid{Dim,T}(dims, origin, spacing, offset) where {Dim,T}
    @assert all(>(0), dims) "dimensions must be positive"
    @assert all(>(0), spacing) "spacing must be positive"
    topology = GridTopology(dims)
    new(origin, spacing, offset, topology)
  end
end

CartesianGrid(dims::Dims{Dim}, origin::Point{Dim,T},
              spacing::NTuple{Dim,T},
              offset::Dims{Dim}=ntuple(i->1, Dim)) where {Dim,T} =
  CartesianGrid{Dim,T}(dims, origin, spacing, offset)

CartesianGrid(dims::Dims{Dim}, origin::NTuple{Dim,T},
              spacing::NTuple{Dim,T},
              offset::Dims{Dim}=ntuple(i->1, Dim)) where {Dim,T} =
  CartesianGrid{Dim,T}(dims, Point(origin), spacing, offset)

function CartesianGrid(start::Point{Dim,T}, finish::Point{Dim,T},
                       spacing::NTuple{Dim,T}) where {Dim,T}
  dims = Tuple(ceil.(Int, (finish - start) ./ spacing))
  origin = start
  offset = ntuple(i->1, Dim)
  CartesianGrid{Dim,T}(dims, origin, spacing, offset)
end

CartesianGrid(start::NTuple{Dim,T}, finish::NTuple{Dim,T},
              spacing::NTuple{Dim,T}) where {Dim,T} =
  CartesianGrid(Point(start), Point(finish), spacing)

function CartesianGrid(start::Point{Dim,T}, finish::Point{Dim,T};
                       dims::Dims{Dim}=ntuple(i->100, Dim)) where {Dim,T}
  origin  = start
  spacing = Tuple((finish - start) ./ dims)
  offset  = ntuple(i->1, Dim)
  CartesianGrid{Dim,T}(dims, origin, spacing, offset)
end

CartesianGrid(start::NTuple{Dim,T}, finish::NTuple{Dim,T};
              dims::Dims{Dim}=ntuple(i->100, Dim)) where {Dim,T} =
  CartesianGrid(Point(start), Point(finish); dims=dims)

function CartesianGrid{T}(dims::Dims{Dim}) where {Dim,T}
  origin  = ntuple(i->zero(T), Dim)
  spacing = ntuple(i->one(T), Dim)
  offset  = ntuple(i->1, Dim)
  CartesianGrid{Dim,T}(dims, origin, spacing, offset)
end

CartesianGrid{T}(dims::Vararg{Int,Dim}) where {Dim,T} = CartesianGrid{T}(dims)

CartesianGrid(dims::Dims{Dim}) where {Dim} = CartesianGrid{Float64}(dims)

CartesianGrid(dims::Vararg{Int,Dim}) where {Dim} = CartesianGrid{Float64}(dims)

Base.size(g::CartesianGrid) = size(g.topology)
Base.minimum(g::CartesianGrid) = Point(coordinates(g.origin) .- (g.offset .- 1) .* g.spacing)
Base.maximum(g::CartesianGrid) = Point(coordinates(g.origin) .+ (size(g.topology) .- g.offset .+ 1) .* g.spacing)
Base.extrema(g::CartesianGrid) = minimum(g), maximum(g)
spacing(g::CartesianGrid) = g.spacing
offset(g::CartesianGrid) = g.offset

==(g1::CartesianGrid, g2::CartesianGrid) =
  g1.topology == g2.topology && g1.spacing  == g2.spacing &&
  Tuple(g1.origin - g2.origin) == (g1.offset .- g2.offset) .* g1.spacing

# -----------------
# DOMAIN INTERFACE
# -----------------

function element(g::CartesianGrid{Dim}, ind::Int) where {Dim}
  dims = size(g.topology)
  I = CartesianIndices(dims)[ind]
  o = coordinates(g.origin)
  s = g.spacing
  i = I.I .- g.offset .+ 1

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

function centroid(g::CartesianGrid{Dim}, ind::Int) where {Dim}
  dims = size(g.topology)
  intcoords = CartesianIndices(dims)[ind]
  neworigin = coordinates(g.origin) .+ g.spacing ./ 2
  Point(ntuple(i -> neworigin[i] + (intcoords[i] - g.offset[i])*g.spacing[i], Dim))
end

Base.eltype(g::CartesianGrid) = typeof(g[1])

# ---------------
# MESH INTERFACE
# ---------------

function vertices(g::CartesianGrid)
  dims = size(g.topology)
  inds = CartesianIndices(dims .+ 1)
  vec([Point(coordinates(g.origin) .+ (ind.I .- g.offset) .* g.spacing) for ind in inds])
end

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
  dims   = size(I)
  offset = g.offset .- first(I).I .+ 1
  CartesianGrid(dims, g.origin, g.spacing, offset)
end

Base.view(g::CartesianGrid{Dim}, I::CartesianIndices{Dim}) where {Dim} = getindex(g, I)

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, g::CartesianGrid{Dim,T}) where {Dim,T}
  dims = join(size(g.topology), "×")
  print(io, "$dims CartesianGrid{$Dim,$T}")
end

function Base.show(io::IO, ::MIME"text/plain", g::CartesianGrid)
  println(io, g)
  println(io, "  minimum: ", minimum(g))
  println(io, "  maximum: ", maximum(g))
  print(  io, "  spacing: ", spacing(g))
end
