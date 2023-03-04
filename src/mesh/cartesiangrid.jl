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
struct CartesianGrid{Dim,T} <: Grid{Dim,T}
  origin::Point{Dim,T}
  spacing::NTuple{Dim,T}
  offset::Dims{Dim}
  topology::GridTopology{Dim}

  function CartesianGrid{Dim,T}(dims, origin, spacing, offset) where {Dim,T}
    @assert all(>(0), dims) "dimensions must be positive"
    @assert all(>(0), spacing) "spacing must be positive"
    topology = GridTopology(dims)
    return new(origin, spacing, offset, topology)
  end
end

function CartesianGrid(
  dims::Dims{Dim},
  origin::Point{Dim,T},
  spacing::NTuple{Dim,T},
  offset::Dims{Dim}=ntuple(i -> 1, Dim),
) where {Dim,T}
  return CartesianGrid{Dim,T}(dims, origin, spacing, offset)
end

function CartesianGrid(
  dims::Dims{Dim},
  origin::NTuple{Dim,T},
  spacing::NTuple{Dim,T},
  offset::Dims{Dim}=ntuple(i -> 1, Dim),
) where {Dim,T}
  return CartesianGrid{Dim,T}(dims, Point(origin), spacing, offset)
end

function CartesianGrid(
  start::Point{Dim,T}, finish::Point{Dim,T}, spacing::NTuple{Dim,T}
) where {Dim,T}
  dims = Tuple(ceil.(Int, (finish - start) ./ spacing))
  origin = start
  offset = ntuple(i -> 1, Dim)
  return CartesianGrid{Dim,T}(dims, origin, spacing, offset)
end

function CartesianGrid(
  start::NTuple{Dim,T}, finish::NTuple{Dim,T}, spacing::NTuple{Dim,T}
) where {Dim,T}
  return CartesianGrid(Point(start), Point(finish), spacing)
end

function CartesianGrid(
  start::Point{Dim,T}, finish::Point{Dim,T}; dims::Dims{Dim}=ntuple(i -> 100, Dim)
) where {Dim,T}
  origin = start
  spacing = Tuple((finish - start) ./ dims)
  offset = ntuple(i -> 1, Dim)
  return CartesianGrid{Dim,T}(dims, origin, spacing, offset)
end

function CartesianGrid(
  start::NTuple{Dim,T}, finish::NTuple{Dim,T}; dims::Dims{Dim}=ntuple(i -> 100, Dim)
) where {Dim,T}
  return CartesianGrid(Point(start), Point(finish); dims=dims)
end

function CartesianGrid{T}(dims::Dims{Dim}) where {Dim,T}
  origin = ntuple(i -> zero(T), Dim)
  spacing = ntuple(i -> one(T), Dim)
  offset = ntuple(i -> 1, Dim)
  return CartesianGrid{Dim,T}(dims, origin, spacing, offset)
end

CartesianGrid{T}(dims::Vararg{Int,Dim}) where {Dim,T} = CartesianGrid{T}(dims)

CartesianGrid(dims::Dims{Dim}) where {Dim} = CartesianGrid{Float64}(dims)

CartesianGrid(dims::Vararg{Int,Dim}) where {Dim} = CartesianGrid{Float64}(dims)

function cart2vert(g::CartesianGrid, ijk::Tuple)
  return Point(coordinates(g.origin) .+ (ijk .- g.offset) .* g.spacing)
end

spacing(g::CartesianGrid) = g.spacing

offset(g::CartesianGrid) = g.offset

function centroid(g::CartesianGrid, ind::Int)
  ijk = elem2cart(topology(g), ind)
  p = cart2vert(g, ijk)
  δ = Vec(spacing(g) ./ 2)
  return p + δ
end

function Base.getindex(g::CartesianGrid{Dim}, I::CartesianIndices{Dim}) where {Dim}
  dims = size(I)
  offset = g.offset .- first(I).I .+ 1
  return CartesianGrid(dims, g.origin, g.spacing, offset)
end

function ==(g1::CartesianGrid, g2::CartesianGrid)
  return g1.topology == g2.topology &&
         g1.spacing == g2.spacing &&
         Tuple(g1.origin - g2.origin) == (g1.offset .- g2.offset) .* g1.spacing
end

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, g::CartesianGrid{Dim,T}) where {Dim,T}
  dims = join(size(g.topology), "×")
  return print(io, "$dims CartesianGrid{$Dim,$T}")
end

function Base.show(io::IO, ::MIME"text/plain", g::CartesianGrid)
  println(io, g)
  println(io, "  minimum: ", minimum(g))
  println(io, "  maximum: ", maximum(g))
  return print(io, "  spacing: ", spacing(g))
end
