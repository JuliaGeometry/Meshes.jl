# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Proj(CRS)

Convert the coordinates of geometry or domain to a given
coordinate reference system `CRS`.

## Examples

```julia
Proj(Polar)
Proj(WebMercator)
Proj(Mercator{WGS84Latest})
```
"""
struct Proj{CRS} <: CoordinateTransform end

Proj(CRS) = Proj{CRS}()

parameters(::Proj{CRS}) where {CRS} = (; CRS)

applycoord(::Proj, v::Vec) = v

applycoord(::Proj{CRS}, p::Point) where {CRS} = Point(convert(CRS, coords(p)))

# --------------
# SPECIAL CASES
# --------------

applycoord(t::Proj, g::RectilinearGrid) = applycoord(t, convert(SimpleMesh, g))

applycoord(t::Proj, g::StructuredGrid) = applycoord(t, convert(SimpleMesh, g))
