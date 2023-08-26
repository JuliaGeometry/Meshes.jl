# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    StdCoords()

Standardize coordinates of all geometries
to the interval `[-0.5, 0.5]`.
"""
struct StdCoords <: GeometricTransform end

isrevertible(::Type{<:StdCoords}) = true

function apply(::StdCoords, g::GeometryOrDomain)
  box = boundingbox(g)
  c, s = center(box), sides(box)
  tr = Translate(coordinates(c)...)
  ts = Stretch(s)
  t = inv(tr) → inv(ts)
  t(g), t
end

revert(::StdCoords, g::GeometryOrDomain, t) = inv(t)(g)

reapply(::StdCoords, g::GeometryOrDomain, t) = t(g)
