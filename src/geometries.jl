# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Geometry{Dim,T}

A geometry embedded in a `Dim`-dimensional space with coordinates of type `T`.
"""
abstract type Geometry{Dim,T} end

Broadcast.broadcastable(g::Geometry) = Ref(g)

"""
    embeddim(geometry)

Return the number of dimensions of the space where the `geometry` is embedded.
"""
embeddim(::Type{<:Geometry{Dim,T}}) where {Dim,T} = Dim
embeddim(g::Geometry) = embeddim(typeof(g))

"""
    paramdim(geometry)

Return the number of parametric dimensions of the `geometry`. For example, a
sphere embedded in 3D has 2 parametric dimensions (polar and azimuthal angles).

See also [`isparametrized`](@ref).
"""
paramdim(g::Geometry) = paramdim(typeof(g))

"""
    coordtype(geometry)

Return the machine type of each coordinate used to describe the `geometry`.
"""
coordtype(::Type{<:Geometry{Dim,T}}) where {Dim,T} = T
coordtype(g::Geometry) = coordtype(typeof(g))

"""
    centroid(geometry)

Return the centroid of the `geometry`.
"""
centroid(g::Geometry) = center(g)

"""
    measure(geometry)

Return the measure of the `geometry`, i.e. the length, area, or volume.
"""
measure(g::Geometry) = sum(measure, discretize(g))

"""
    perimeter(geometry)

Return the perimeter of the `geometry`, i.e. the measure of its boundary.
"""
perimeter(g::Geometry) = measure(boundary(g))

"""
    extrema(geometry)

Return the top left and bottom right corners of the
bounding box of the `geometry`.
"""
Base.extrema(g::Geometry) = extrema(boundingbox(g))

"""
    boundary(geometry)

Return the boundary of the `geometry`.
"""
function boundary end

"""
    g₁ ⊆ g₂

Tells whether or not the geometry `g₁` is a subset of geometry `g₂`.
"""
Base.issubset(g₁::Geometry, g₂::Geometry) = all(p ∈ g₂ for p in vertices(g₁))

"""
    isconvex(geometry)

Tells whether or not the `geometry` is convex.
"""
isconvex(g::Geometry) = isconvex(typeof(g))

"""
    issimplex(geometry)

Tells whether or not the `geometry` is simplex.
"""
issimplex(::Type{<:Geometry}) = false
issimplex(g::Geometry) = issimplex(typeof(g))

"""
    isperiodic(geometry)

Tells whether or not the `geometry` is periodic
along each parametric dimension.
"""
isperiodic(g::Geometry) = isperiodic(typeof(g))

"""
    isparametrized(geometry)

Tells whether or not the `geometry` is parametrized,
i.e. can be called as `geometry(u₁, u₂, ..., uₙ)` with
local coordinates `(u₁, u₂, ..., uₙ) ∈ [0,1]ⁿ` where
`n` is the parametric dimension.

See also [`paramdim`](@ref).
"""
isparametrized(::Type{<:Geometry}) = false
isparametrized(g::Geometry) = isparametrized(typeof(g))

# ----------------
# IMPLEMENTATIONS
# ----------------

include("primitives.jl")
include("polytopes.jl")
include("multigeoms.jl")
