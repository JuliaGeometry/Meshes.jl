# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# alias for dispatch purposes
const QuasiCartesianGrid{M<:𝔼,C<:Union{Cartesian,Projected}} = RegularGrid{M,C}

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

`CartesianGrid` is an alias to [`RegularGrid`](@ref) with `Cartesian` CRS.

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

See also [`RegularGrid`](@ref).
"""
const CartesianGrid{M<:𝔼,C<:Cartesian} = RegularGrid{M,C}

CartesianGrid(
  origin::Point{𝔼{Dim}},
  spacing::NTuple{Dim,Number},
  offset::Dims{Dim},
  topology::GridTopology{Dim}
) where {Dim} = RegularGrid(_cartpoint(origin), spacing, offset, topology)

CartesianGrid(
  dims::Dims{Dim},
  origin::Point{𝔼{Dim}},
  spacing::NTuple{Dim,Number},
  offset::Dims{Dim}=ntuple(i -> 1, Dim)
) where {Dim} = RegularGrid(dims, _cartpoint(origin), spacing, offset)

CartesianGrid(
  dims::Dims{Dim},
  origin::NTuple{Dim,Number},
  spacing::NTuple{Dim,Number},
  offset::Dims{Dim}=ntuple(i -> 1, Dim)
) where {Dim} = CartesianGrid(dims, Point(origin), spacing, offset)

CartesianGrid(start::Point{𝔼{Dim}}, finish::Point{𝔼{Dim}}, spacing::NTuple{Dim,Number}) where {Dim} =
  RegularGrid(_cartpoint(start), _cartpoint(finish), spacing)

CartesianGrid(start::NTuple{Dim,Number}, finish::NTuple{Dim,Number}, spacing::NTuple{Dim,Number}) where {Dim} =
  CartesianGrid(Point(start), Point(finish), spacing)

CartesianGrid(start::Point{𝔼{Dim}}, finish::Point{𝔼{Dim}}; dims::Dims{Dim}=ntuple(i -> 100, Dim)) where {Dim} =
  RegularGrid(_cartpoint(start), _cartpoint(finish); dims)

CartesianGrid(
  start::NTuple{Dim,Number},
  finish::NTuple{Dim,Number};
  dims::Dims{Dim}=ntuple(i -> 100, Dim)
) where {Dim} = CartesianGrid(Point(start), Point(finish); dims)

function CartesianGrid(dims::Dims{Dim}) where {Dim}
  origin = ntuple(i -> 0.0, Dim)
  spacing = ntuple(i -> 1.0, Dim)
  offset = ntuple(i -> 1, Dim)
  CartesianGrid(dims, origin, spacing, offset)
end

CartesianGrid(dims::Int...) = CartesianGrid(dims)

# -----------------
# HELPER FUNCTIONS
# -----------------

_cartpoint(p) = Point(convert(Cartesian, coords(p)))
