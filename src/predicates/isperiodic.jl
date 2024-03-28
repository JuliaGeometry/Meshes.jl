# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    isperiodic(geometry)

Tells whether or not the `geometry` is periodic
along each parametric dimension.
"""
isperiodic(g::Geometry) = isperiodic(typeof(g))

isperiodic(::Type{<:Segment}) = (false,)

isperiodic(::Type{<:Ray}) = (false,)

isperiodic(::Type{<:Line}) = (false,)

isperiodic(b::BezierCurve) = (first(controls(b)) == last(controls(b)),)

isperiodic(::Type{<:Plane}) = (false, false)

isperiodic(::Type{<:Box{Dim}}) where {Dim} = ntuple(i -> false, Dim)

isperiodic(::Type{<:Ball{Dim}}) where {Dim} = ntuple(i -> i != 1, Dim)

isperiodic(::Type{<:Sphere{Dim}}) where {Dim} = ntuple(i -> true, Dim - 1)

isperiodic(::Type{<:Disk}) = (false, true)

isperiodic(::Type{<:Circle}) = (true,)

isperiodic(::Type{<:CylinderSurface}) = (true, false)

isperiodic(::Type{<:ConeSurface}) = (true, false)

isperiodic(::Type{<:FrustumSurface}) = (true, false)

isperiodic(::Type{<:ParaboloidSurface}) = (false, true)

isperiodic(::Type{<:Torus}) = (true, true)

isperiodic(c::Type{<:Chain}) = (isclosed(c),)

isperiodic(::Type{<:Quadrangle}) = (false, false)

isperiodic(::Type{<:Hexahedron}) = (false, false, false)

"""
    isperiodic(grid)

Tells whether or not the `grid` is periodic
along each parametric dimension.
"""
isperiodic(g::Grid) = isperiodic(topology(g))
