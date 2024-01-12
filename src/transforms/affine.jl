# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Affine(rot, offsets...)
    Affine((u, v), offsets...)
    Affine(θ, offsets...)

Rotate geometry or mesh with rotation `rot` and translate
their coordinates by given offsets `o₁, o₂, ...`.
This transform is equivalet to `Rotate(rot) → Translate(offsets...)`.

The rotation argument can also be a tuple of vectors `(u, v)` or a
rotation angle `θ`. See [`Rotate`](@ref) documentation for more details.

# Examples

```julia
Affine(AngleAxis(0.2, 1.0, 0.0, 0.0), -2, 2, 2)
Affine((Vec(1, 0, 0), Vec(1, 1, 1)), 2, -2, 2)
Affine(π / 2, 2, -2)
```

See also [`Rotate`](@ref), [`Translate`](@ref).
"""
struct Affine{Dim,T,R<:Rotation{Dim,T}} <: CoordinateTransform
  rot::R
  offsets::NTuple{Dim,T}
end

Affine(rot::R, offsets::Tuple) where {Dim,T,R<:Rotation{Dim,T}} = Affine{Dim,T,R}(rot, offsets)

Affine(rot::Rotation, offsets...) = Affine(rot, offsets)

Affine((u, v)::NTuple{2,Vec}, offsets...) = Affine(rotation_between(u, v), offsets)
Affine((u, v)::NTuple{2,Tuple}, offsets...) = Affine((Vec(u), Vec(v)), offsets...)

Affine(θ, offsets...) = Affine(Angle2d(θ), offsets)

isrevertible(::Type{<:Affine}) = false

isinvertible(::Type{<:Affine}) = false

applycoord(t::Affine, v::Vec) = t.rot * v

applycoord(t::Affine, p::Point) = Point(t.rot * coordinates(p)) + Vec(t.offsets)

# --------------
# SPECIAL CASES
# --------------

applycoord(t::Affine, b::Box{2}) = applycoord(t, convert(Quadrangle, b))

applycoord(t::Affine, b::Box{3}) = applycoord(t, convert(Hexahedron, b))
