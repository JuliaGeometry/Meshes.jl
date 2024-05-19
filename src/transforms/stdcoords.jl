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

function apply(t::StdCoords, g::GeometryOrDomain)
  p = _stdcoords(t, g)
  n, c = apply(p, g)
  n, (p, c)
end

revert(t::StdCoords, g::GeometryOrDomain, c) = revert(c[1], g, c[2])

reapply(t::StdCoords, g::GeometryOrDomain, c) = reapply(c[1], g, c[2])

function _stdcoords(t, g)
  b = boundingbox(g)
  t = Translate(coordinates(center(b))...)
  s = Scale(ustrip.(sides(b)))
  inverse(t) → inverse(s)
end
