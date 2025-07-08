# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    isparametrized(object)

Tells whether or not the geometric `object` is parametrized, i.e.
can be called as `object(u₁, u₂, ..., uₙ)` with local coordinates
`(u₁, u₂, ..., uₙ) ∈ [0,1]ⁿ` where `n` is the parametric dimension.

See also [`paramdim`](@ref).
"""
isparametrized(g) = isparametrized(typeof(g))

isparametrized(::Type{<:Geometry}) = false

isparametrized(::Type{<:Segment}) = true

isparametrized(::Type{<:Ray}) = true

isparametrized(::Type{<:Line}) = true

isparametrized(::Type{<:Plane}) = true

isparametrized(::Type{<:BezierCurve}) = true

isparametrized(::Type{<:ParametrizedCurve}) = true

isparametrized(::Type{<:Box{<:𝔼}}) = true

isparametrized(::Type{<:Ball{<:𝔼}}) = true

isparametrized(::Type{<:Sphere{<:𝔼}}) = true

isparametrized(::Type{<:Ellipsoid}) = true

isparametrized(::Type{<:Disk}) = true

isparametrized(::Type{<:Circle}) = true

isparametrized(::Type{<:Cylinder}) = true

isparametrized(::Type{<:CylinderSurface}) = false

isparametrized(::Type{<:CylinderWall}) = true

isparametrized(::Type{<:Cone}) = true

isparametrized(::Type{<:ConeSurface}) = true

isparametrized(::Type{<:FrustumSurface}) = true

isparametrized(::Type{<:ParaboloidSurface}) = true

isparametrized(::Type{<:Torus}) = true

isparametrized(::Type{<:Triangle}) = true

isparametrized(::Type{<:Quadrangle}) = true

isparametrized(::Type{<:Tetrahedron}) = true

isparametrized(::Type{<:Hexahedron}) = true

isparametrized(::Type{<:TransformedGeometry{M,C,G}}) where {M,C,G} = isparametrized(G)

isparametrized(::Type{<:Domain}) = false
