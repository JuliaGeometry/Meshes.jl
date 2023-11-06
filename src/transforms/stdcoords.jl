# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    StdCoords()

Standardize coordinates of all geometries
to the interval `[-0.5, 0.5]`.

## Examples

```julia
julia> CartesianGrid(10, 10) |> StdCoords()
10×10 CartesianGrid{2,Float64}
  minimum: Point(-0.5, -0.5)
  maximum: Point(0.5, 0.5)
  spacing: (0.1, 0.1)
```
"""
struct StdCoords <: GeometricTransform end

isrevertible(::Type{<:StdCoords}) = true

function apply(::StdCoords, g::GeometryOrDomain)
  box = boundingbox(g)
  c, s = center(box), sides(box)
  tr = Translate(coordinates(c)...)
  ts = Stretch(s)
  t = inverse(tr) → inverse(ts)
  t(g), t
end

revert(::StdCoords, g::GeometryOrDomain, t) = inverse(t)(g)

reapply(::StdCoords, g::GeometryOrDomain, t) = t(g)
