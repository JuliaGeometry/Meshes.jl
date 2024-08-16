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

# convert the CRS and use the default manifold
applycoord(::Proj{CRS}, p::Point) where {CRS} = Point(convert(CRS, coords(p)))

# convert the CRS and preserve the manifold
applycoord(::Proj{CRS}, p::Point{M}) where {CRS<:Basic,M<:🌐} = Point{M}(convert(CRS, coords(p)))

# fix ambiguities
applycoord(::Proj{CRS}, p::Point{<:🌐}) where {CRS<:Projected} = Point(convert(CRS, coords(p)))
applycoord(::Proj{CRS}, p::Point{<:𝔼}) where {CRS<:Geographic} = Point(convert(CRS, coords(p)))

applycoord(::Proj, v::Vec) = v

# --------------
# SPECIAL CASES
# --------------

applycoord(t::Proj{<:Projected}, g::Primitive{<:🌐}) = TransformedGeometry(g, t)

applycoord(t::Proj{<:Geographic}, g::Primitive{<:𝔼}) = TransformedGeometry(g, t)

applycoord(t::Proj, g::RectilinearGrid) = applycoord(t, convert(SimpleMesh, g))

applycoord(t::Proj, g::StructuredGrid) = applycoord(t, convert(SimpleMesh, g))

# -----------
# IO METHODS
# -----------

Base.show(io::IO, ::Proj{CRS}) where {CRS} = print(io, "Proj(CRS: $CRS)")

function Base.show(io::IO, ::MIME"text/plain", t::Proj{CRS}) where {CRS}
  summary(io, t)
  println(io)
  print(io, "└─ CRS: $CRS")
end
