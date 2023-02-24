# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Translate(o₁, o₂, ...)

Translate coordinates of geometry or mesh by
given offsets `o₁, o₂, ...`.
"""
struct Translate{Dim,T} <: StatelessGeometricTransform
  offsets::NTuple{Dim,T}
end

Translate(offsets...) = Translate(offsets)

isrevertible(::Type{<:Translate}) = true

preprocess(transform::Translate, object) = transform.offsets

function applypoint(::Translate, points, prep)
  o = prep
  newpoints = [Point(coordinates(p) .+ o) for p in points]
  newpoints, prep
end

function revertpoint(::Translate, newpoints, cache)
  o = cache
  [Point(coordinates(p) .- o) for p in newpoints]
end