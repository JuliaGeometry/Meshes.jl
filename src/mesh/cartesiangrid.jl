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
julia> CartesianGrid(100, 100, 50)
```

Create a 2D grid with 100 x 100 quadrangles and origin at (10.0, 20.0):

```julia
julia> CartesianGrid((100, 100), (10.0, 20.0), (1.0, 1.0))
```

Create a 1D grid from -1 to 1 with 100 segments:

```julia
julia> CartesianGrid((-1.0,), (1.0,), dims=(100,))
```
"""
struct CartesianGrid{Dim,T} <: Grid{Dim,T}
  origin::Point{Dim,T}
  spacing::NTuple{Dim,T}
  offset::Dims{Dim}
  topology::GridTopology{Dim}
end

function CartesianGrid(
  dims::Dims{Dim},
  origin::Point{Dim,T},
  spacing::NTuple{Dim,T},
  offset::Dims{Dim}=ntuple(i -> 1, Dim)
) where {Dim,T}
  @assert all(>(0), dims) "dimensions must be positive"
  @assert all(>(zero(T)), spacing) "spacing must be positive"
  CartesianGrid{Dim,T}(origin, spacing, offset, GridTopology(dims))
end

CartesianGrid(
  dims::Dims{Dim},
  origin::NTuple{Dim,T},
  spacing::NTuple{Dim,T},
  offset::Dims{Dim}=ntuple(i -> 1, Dim)
) where {Dim,T} = CartesianGrid(dims, Point(origin), spacing, offset)

function CartesianGrid(start::Point{Dim,T}, finish::Point{Dim,T}, spacing::NTuple{Dim,T}) where {Dim,T}
  dims = Tuple(ceil.(Int, (finish - start) ./ spacing))
  origin = start
  offset = ntuple(i -> 1, Dim)
  CartesianGrid(dims, origin, spacing, offset)
end

CartesianGrid(start::NTuple{Dim,T}, finish::NTuple{Dim,T}, spacing::NTuple{Dim,T}) where {Dim,T} =
  CartesianGrid(Point(start), Point(finish), spacing)

function CartesianGrid(start::Point{Dim,T}, finish::Point{Dim,T}; dims::Dims{Dim}=ntuple(i -> 100, Dim)) where {Dim,T}
  origin = start
  spacing = Tuple((finish - start) ./ dims)
  offset = ntuple(i -> 1, Dim)
  CartesianGrid(dims, origin, spacing, offset)
end

CartesianGrid(start::NTuple{Dim,T}, finish::NTuple{Dim,T}; dims::Dims{Dim}=ntuple(i -> 100, Dim)) where {Dim,T} =
  CartesianGrid(Point(start), Point(finish); dims=dims)

function CartesianGrid{T}(dims::Dims{Dim}) where {Dim,T}
  origin = ntuple(i -> zero(T), Dim)
  spacing = ntuple(i -> oneunit(T), Dim)
  offset = ntuple(i -> 1, Dim)
  CartesianGrid(dims, origin, spacing, offset)
end

CartesianGrid{T}(dims::Vararg{Int,Dim}) where {Dim,T} = CartesianGrid{T}(dims)

CartesianGrid(dims::Dims{Dim}) where {Dim} = CartesianGrid{Float64}(dims)

CartesianGrid(dims::Vararg{Int,Dim}) where {Dim} = CartesianGrid{Float64}(dims)

vertex(g::CartesianGrid{Dim}, ijk::Dims{Dim}) where {Dim} =
  Point(coordinates(g.origin) .+ (ijk .- g.offset) .* g.spacing)

spacing(g::CartesianGrid) = g.spacing

offset(g::CartesianGrid) = g.offset

function xyz(g::CartesianGrid{Dim}) where {Dim}
  dims = size(g)
  spac = spacing(g)
  orig = coordinates(minimum(g))
  ntuple(Dim) do i
    o, s, d = orig[i], spac[i], dims[i]
    range(start=o, step=s, length=(d + 1))
  end
end

XYZ(g::CartesianGrid) = XYZ(xyz(g))

function centroid(g::CartesianGrid, ind::Int)
  ijk = elem2cart(topology(g), ind)
  p = vertex(g, ijk)
  δ = Vec(spacing(g) ./ 2)
  p + δ
end

function Base.getindex(g::CartesianGrid{Dim}, I::CartesianIndices{Dim}) where {Dim}
  @boundscheck _checkbounds(g, I)
  dims = size(I)
  offset = g.offset .- Tuple(first(I)) .+ 1
  CartesianGrid(dims, g.origin, g.spacing, offset)
end

==(g1::CartesianGrid, g2::CartesianGrid) =
  g1.topology == g2.topology &&
  g1.spacing == g2.spacing &&
  Tuple(g1.origin - g2.origin) == (g1.offset .- g2.offset) .* g1.spacing

# -----------
# IO METHODS
# -----------

function Base.summary(io::IO, g::CartesianGrid{Dim,T}) where {Dim,T}
  dims = join(size(g.topology), "×")
  print(io, "$dims CartesianGrid{$Dim,$T}")
end

Base.show(io::IO, g::CartesianGrid) = summary(io, g)

function Base.show(io::IO, ::MIME"text/plain", g::CartesianGrid)
  println(io, g)
  println(io, "  minimum: ", minimum(g))
  println(io, "  maximum: ", maximum(g))
  print(io, "  spacing: ", spacing(g))
end
