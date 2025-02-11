# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Proj(CRS)
    Proj(code)

Convert the coordinates of geometry or domain to a given
coordinate reference system `CRS` or EPSG/ESRI `code`.

## Examples

```julia
Proj(Polar)
Proj(WebMercator)
Proj(Mercator{WGS84Latest})
Proj(EPSG{3395})
Proj(ESRI{54017})
```
"""
struct Proj{CRS} <: CoordinateTransform end

Proj(CRS) = Proj{CRS}()

Proj(code::Type{<:EPSG}) = Proj{CoordRefSystems.get(code)}()

Proj(code::Type{<:ESRI}) = Proj{CoordRefSystems.get(code)}()

parameters(::Proj{CRS}) where {CRS} = (; CRS)

applycoord(::Proj, v::Vec) = v

applycoord(::Proj{CRS}, p::Point) where {CRS} = Point(convert(CRS, coords(p)))

# --------------
# SPECIAL CASES
# --------------

applycoord(t::Proj, g::RectilinearGrid) = applycoord(t, convert(SimpleMesh, g))

applycoord(t::Proj, g::StructuredGrid) = applycoord(t, convert(SimpleMesh, g))
