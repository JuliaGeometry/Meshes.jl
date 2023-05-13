# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    StdCoords()

Transform points to lie within box of
unitary sides centered at the origin.
"""
struct StdCoords <: GeometricTransform end

isrevertible(::Type{<:StdCoords}) = true

function preprocess(::StdCoords, object)
  pset = PointSet(_points(object))
  bbox = boundingbox(pset)
  center(bbox), sides(bbox)
end

function applypoint(::StdCoords, points, prep)
  pₒ, s = prep
  newpoints = [Point((p - pₒ) ./ s) for p in points]
  newpoints, prep
end

function revertpoint(::StdCoords, newpoints, pcache)
  pₒ, s = pcache
  cₒ = coordinates(pₒ)
  [Point(s .* coordinates(p) + cₒ) for p in newpoints]
end

function reapplypoint(::StdCoords, points, pcache)
  pₒ, s = pcache
  [Point((p - pₒ) ./ s) for p in points]
end
