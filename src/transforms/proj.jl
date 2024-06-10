# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Proj(CRS)

Projects the coordinates of geometry or domain into `CRS`.

## Examples

```julia
Proj(Polar)
Proj(Mercator{WGS84Latest})
```
"""
struct Proj{CRS} <: CoordinateTransform end

Proj(CRS) = Proj{CRS}()

parameters(::Proj{CRS}) where {CRS} = (; CRS)

applycoord(::Proj, v::Vec) = v

applycoord(::Proj{CRS}, p::Point) where {CRS} = Point(convert(CRS, coords(p)))
