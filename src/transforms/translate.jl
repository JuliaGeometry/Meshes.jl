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
  Translate(offsets::NTuple{Dim,ℒ}) where {Dim,ℒ<:Len} = new{Dim,float(ℒ)}(offsets)
end

Translate(offsets::NTuple{Dim,Len}) where {Dim} = Translate(promote(offsets...))

Translate(offsets::Tuple) = Translate(addunit.(offsets, u"m"))

Translate(offsets...) = Translate(offsets)

parameters(t::Translate) = (; offsets=t.offsets)

isaffine(::Type{<:Translate}) = true

isrevertible(::Type{<:Translate}) = true

isinvertible(::Type{<:Translate}) = true

inverse(t::Translate) = Translate(-1 .* t.offsets)

applycoord(t::Translate, p::Point) = p + Vec(t.offsets)

applycoord(::Translate, v::Vec) = v

# --------------
# SPECIAL CASES
# --------------

apply(t::Translate, g::RectilinearGrid) =
  RectilinearGrid{manifold(g),crs(g)}(ntuple(i -> xyz(g)[i] .+ t.offsets[i], paramdim(g))), nothing

revert(t::Translate, g::RectilinearGrid, c) = first(apply(inverse(t), g))

apply(t::Translate, g::StructuredGrid) =
  StructuredGrid{manifold(g),crs(g)}(ntuple(i -> XYZ(g)[i] .+ t.offsets[i], paramdim(g))), nothing

revert(t::Translate, g::StructuredGrid, c) = first(apply(inverse(t), g))
