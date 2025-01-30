# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    CartesianGrid(args...; kwargs...)

A Cartesian grid is a [`RegularGrid`](@ref) where all arguments
are forced to have `Cartesian` coordinates. Please check the
docstring of [`RegularGrid`](@ref) for more information on
possible `args` and `kwargs`.

See also [`RegularGrid`](@ref).
"""
const CartesianGrid{M<:𝔼,C<:Cartesian} = RegularGrid{M,C}

CartesianGrid(args...; kwargs...) = RegularGrid(_cartesian.(args)...; kwargs...)

# enforce Cartesian coordinates for all Point arguments
_cartesian(p::Point) = Point(convert(Cartesian, coords(p)))

# forward all other arguments without change
_cartesian(o) = o
