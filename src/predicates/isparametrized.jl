# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    isparametrized(geometry)

Tells whether or not the `geometry` is parametrized,
i.e. can be called as `geometry(u₁, u₂, ..., uₙ)` with
local coordinates `(u₁, u₂, ..., uₙ) ∈ [0,1]ⁿ` where
`n` is the parametric dimension.

See also [`paramdim`](@ref).
"""
function isparametrized end

isparametrized(g::Geometry) = isparametrized(typeof(g))

isparametrized(::Type{<:Geometry}) = false

isparametrized(::Type{<:Segment}) = true

isparametrized(::Type{<:Ray}) = true

isparametrized(::Type{<:Line}) = true

isparametrized(::Type{<:Plane}) = true

isparametrized(::Type{<:BezierCurve}) = true

isparametrized(::Type{<:Box}) = true

isparametrized(::Type{<:Ball}) = true

isparametrized(::Type{<:Sphere}) = true

isparametrized(::Type{<:Ellipsoid}) = true

isparametrized(::Type{<:Disk}) = true

isparametrized(::Type{<:Circle}) = true

isparametrized(::Type{<:Cylinder}) = true

isparametrized(::Type{<:CylinderSurface}) = true

isparametrized(::Type{<:ConeSurface}) = true

isparametrized(::Type{<:FrustumSurface}) = true

isparametrized(::Type{<:ParaboloidSurface}) = true

isparametrized(::Type{<:Torus}) = true

isparametrized(::Type{<:Triangle}) = true

isparametrized(::Type{<:Quadrangle}) = true

isparametrized(::Type{<:Tetrahedron}) = true

isparametrized(::Type{<:Hexahedron}) = true

isparametrized(d::Domain) = isparametrized(typeof(d))

isparametrized(::Type{<:Domain}) = false
