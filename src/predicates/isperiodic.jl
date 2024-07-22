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

isperiodic(B::Type{<:Box}) = ntuple(i -> false, embeddim(B))

isperiodic(B::Type{<:Ball}) = ntuple(i -> i == embeddim(B), embeddim(B))

isperiodic(S::Type{<:Sphere}) = ntuple(i -> i == embeddim(S) - 1, embeddim(S) - 1)

isperiodic(::Type{<:Ellipsoid}) = (false, true)

isperiodic(::Type{<:Disk}) = (false, true)

isperiodic(::Type{<:Circle}) = (true,)

isperiodic(::Type{<:Cylinder}) = (false, true, false)

isperiodic(::Type{<:CylinderSurface}) = (true, false)

isperiodic(::Type{<:ConeSurface}) = (true, false)

isperiodic(::Type{<:FrustumSurface}) = (true, false)

isperiodic(::Type{<:ParaboloidSurface}) = (false, true)

isperiodic(::Type{<:Torus}) = (true, true)

isperiodic(::Type{<:Rope}) = (false,)

isperiodic(::Type{<:Ring}) = (true,)

isperiodic(::Type{<:Quadrangle}) = (false, false)

isperiodic(::Type{<:Hexahedron}) = (false, false, false)

"""
    isperiodic(grid)

Tells whether or not the `grid` is periodic
along each parametric dimension.
"""
isperiodic(g::Grid) = isperiodic(topology(g))
