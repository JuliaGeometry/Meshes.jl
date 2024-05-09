# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Translate(o₁, o₂, ...)

Translate coordinates of geometry or mesh by
given offsets `o₁, o₂, ...`.
"""
struct Translate{Dim,ℒ<:Len} <: CoordinateTransform
  offsets::NTuple{Dim,ℒ}
end

# TODO: convert `offsets` values to float?
Translate(offsets::Tuple) = Translate(addunit.(offsets, u"m"))

Translate(offsets...) = Translate(offsets)

parameters(t::Translate) = (; offsets=t.offsets)

isaffine(::Type{<:Translate}) = true

isrevertible(::Type{<:Translate}) = true

isinvertible(::Type{<:Translate}) = true

inverse(t::Translate) = Translate(-1 .* t.offsets)

applycoord(t::Translate, v::Vec) = v

# TODO: should offsets have units?
applycoord(t::Translate, p::Point) = p + Vec(t.offsets)

# ----------------
# SPECIALIZATIONS
# ----------------

apply(t::Translate{Dim}, g::RectilinearGrid{Dim}) where {Dim} =
  RectilinearGrid(ntuple(i -> xyz(g)[i] .+ t.offsets[i], Dim)), nothing

revert(t::Translate, g::RectilinearGrid, c) = first(apply(inverse(t), g))

apply(t::Translate{Dim}, g::StructuredGrid{Dim}) where {Dim} =
  StructuredGrid(ntuple(i -> XYZ(g)[i] .+ t.offsets[i], Dim)), nothing

revert(t::Translate, g::StructuredGrid, c) = first(apply(inverse(t), g))
