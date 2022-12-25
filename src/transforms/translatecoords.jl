# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    TranslateCoords(o₁, o₂, ...)

Translate coordinates of geometry or mesh by
given offsets `o₁, o₂, ...`.
"""
struct TranslateCoords{Dim,T} <: StatelessGeometricTransform
  offsets::NTuple{Dim,T}
end

TranslateCoords(offsets::NTuple{Dim,T}) where {Dim,T} =
  TranslateCoords{Dim,T}(offsets)

TranslateCoords(offsets...) = TranslateCoords(offsets)

isrevertible(::Type{<:TranslateCoords}) = true

preprocess(transform::TranslateCoords, object) = transform.offsets

function applypoint(::TranslateCoords, points, prep)
  o = prep
  newpoints = [Point(coordinates(p) .+ o) for p in points]
  newpoints, prep
end

function revertpoint(::TranslateCoords, newpoints, cache)
  o = cache
  [Point(coordinates(p) .- o) for p in newpoints]
end