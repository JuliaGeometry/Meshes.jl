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

Base.inv(t::Translate) = Translate(-1 .* t.offsets)

applycoord(t::Translate, v::Vec) = v + Vec(t.offsets)

# --------------
# SPECIAL CASES
# --------------

function applycoord(t::Translate, g::CartesianGrid)
  dims = size(g)
  orig = applycoord(t, minimum(g))
  spac = spacing(g)
  offs = offset(g)
  CartesianGrid(dims, orig, spac, offs)
end
