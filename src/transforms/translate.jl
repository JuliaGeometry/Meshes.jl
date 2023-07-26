# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Translate(o₁, o₂, ...)

Translate coordinates of geometry or mesh by
given offsets `o₁, o₂, ...`.
"""
struct Translate{Dim,T} <: CoordinateTransform
  offsets::NTuple{Dim,T}
end

Translate(offsets...) = Translate(offsets)

Base.inv(t::Translate) = Translate(-1 .* t.offsets)

isrevertible(::Type{<:Translate}) = true

_apply(t::Translate, v::Vec) = v + Vec(t.offsets)
