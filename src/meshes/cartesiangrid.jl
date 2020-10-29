# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CartesianGrid(dims, origin, spacing)

A Cartesian grid with dimensions `dims`, lower left corner at `origin`
and cell spacing `spacing`. The three arguments must have the same length.

    CartesianGrid(start, finish, dims=dims)

Alternatively, construct a Cartesian grid from a `start` point (lower left)
to a `finish` point (upper right).

    CartesianGrid{T}(dims)
    CartesianGrid{T}(dim1, dim2, ...)

Finally, a Cartesian grid can be constructed by only passing the dimensions
`dims` as a tuple, or by passing each dimension `dim1`, `dim2`, ... separately.
In this case, the origin and spacing default to (0,0,...) and (1,1,...).

## Examples

Create a 3D grid with 100x100x50 locations:

```julia
julia> CartesianGrid(100,100,50)
```

Create a 2D grid with 100x100 locations and origin at (10.,20.) units:

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
  origin::Point{Dim,T}
  spacing::SVector{Dim,T}

  function CartesianGrid{Dim,T}(dims, origin, spacing) where {Dim,T}
    @assert all(dims .> 0) "dimensions must be positive"
    @assert all(spacing .> 0) "spacing must be positive"
    new(dims, origin, spacing)
  end
end

CartesianGrid(dims::Dims{Dim}, origin::Point{Dim,T}, spacing::SVector{Dim,T}) where {Dim,T} =
  CartesianGrid{Dim,T}(dims, origin, spacing)

CartesianGrid(dims::Dims{Dim}, origin::NTuple{Dim,T}, spacing::NTuple{Dim,T}) where {Dim,T} =
  CartesianGrid{Dim,T}(dims, Point(origin), SVector(spacing))

CartesianGrid(start::Point{Dim,T}, finish::Point{Dim,T};
            dims::Dims{Dim}=ntuple(i->100, Dim)) where {Dim,T} =
  CartesianGrid{Dim,T}(dims, start, (finish - start) ./ (dims .- 1))

CartesianGrid(start::NTuple{Dim,T}, finish::NTuple{Dim,T};
            dims::Dims{Dim}=ntuple(i->100, Dim)) where {Dim,T} =
  CartesianGrid(Point(start), Point(finish); dims=dims)

CartesianGrid{T}(dims::Dims{Dim}) where {Dim,T} =
  CartesianGrid{Dim,T}(dims, ntuple(i->zero(T), Dim), ntuple(i->one(T), Dim))

CartesianGrid{T}(dims::Vararg{Int,Dim}) where {Dim,T} = CartesianGrid{T}(dims)

CartesianGrid(dims::Dims{Dim}) where {Dim} = CartesianGrid{Float64}(dims)

CartesianGrid(dims::Vararg{Int,Dim}) where {Dim} = CartesianGrid{Float64}(dims)

Base.size(g::CartesianGrid) = g.dims
Base.minimum(g::CartesianGrid) = g.origin
Base.maximum(g::CartesianGrid) = g.origin + (g.dims .- 1) .* g.spacing
spacing(g::CartesianGrid) = g.spacing

function Base.show(io::IO, g::CartesianGrid{Dim,T}) where {Dim,T}
  dims = join(g.dims, "×")
  print(io, "$dims CartesianGrid{$Dim,$T}")
end

function Base.show(io::IO, ::MIME"text/plain", g::CartesianGrid{Dim,T}) where {Dim,T}
  println(io, g)
  println(io, "  minimum: ", minimum(g))
  println(io, "  maximum: ", maximum(g))
  print(  io, "  spacing: ", spacing(g))
end
