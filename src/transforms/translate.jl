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

isrevertible(::Type{<:Translate}) = true

isinvertible(::Type{<:Translate}) = true

inverse(t::Translate) = Translate(-1 .* t.offsets)

applycoord(t::Translate, v::Vec) = v

applycoord(t::Translate, p::Point) = p + Vec(t.offsets)
