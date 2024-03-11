# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    issimplex(geometry)

Tells whether or not the `geometry` is a simplex.
"""
issimplex(g::Geometry) = issimplex(typeof(g))

issimplex(::Type{<:Geometry}) = false

issimplex(::Type{<:Triangle}) = true

issimplex(::Type{<:Tetrahedron}) = true

"""
    issimplex(connectivity)

Tells whether or not the `connectivity` is a simplex.
"""
issimplex(c::Connectivity) = issimplex(typeof(c))

issimplex(::Type{Connectivity{PL,N}}) where {PL,N} = issimplex(PL)
