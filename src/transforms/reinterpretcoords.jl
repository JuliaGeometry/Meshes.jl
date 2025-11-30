# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ReinterpretCoords(CRS₁, CRS₂)

Reinterpret the coordinate reference system `CRS₁`
as a different coordinate reference system `CRS₂`
with the same number of coordinates.

    ReinterpretCoords(code₁, code₂)

Alternatively, specify coordinate reference systems
with EPSG/ESRI `code₁` and `code₂`.

## Examples

```julia
ReinterpretCoords(Cartesian, LatLon)
ReinterpretCoords(LatLon, Mercator)
```
"""
struct ReinterpretCoords{CRS₁,CRS₂} <: CoordinateTransform end

ReinterpretCoords(CRS₁, CRS₂) = ReinterpretCoords{CRS₁,CRS₂}()

ReinterpretCoords(code₁::Type{<:CoordRefSystems.CRSCode}, code₂::Type{<:CoordRefSystems.CRSCode}) =
  ReinterpretCoords(CoordRefSystems.get(code₁), CoordRefSystems.get(code₂))

parameters(::ReinterpretCoords{CRS₁,CRS₂}) where {CRS₁,CRS₂} = (; CRS₁, CRS₂)

isrevertible(::Type{<:ReinterpretCoords}) = true

isinvertible(::Type{<:ReinterpretCoords}) = true

inverse(::ReinterpretCoords{CRS₁,CRS₂}) where {CRS₁,CRS₂} = ReinterpretCoords(CRS₂, CRS₁)

function applycoord(::ReinterpretCoords{CRS₁,CRS₂}, p::Point{M,<:CRS₁}) where {M,CRS₁,CRS₂}
  rawcoords = CoordRefSystems.raw(coords(p))
  newcoords = CoordRefSystems.reconstruct(CRS₂, rawcoords)
  Point(newcoords)
end

function applycoord(::ReinterpretCoords{CRS₁,CRS₂}, p::Point{M,CRS}) where {M,CRS₁,CRS₂,CRS}
  throw(ArgumentError("""
    attempting to reinterpret $CRS₁ as $CRS₂ with geometry that has a different CRS:

    $CRS
    """))
end

applycoord(::ReinterpretCoords, v::Vec) = v

# --------------
# SPECIAL CASES
# --------------

applycoord(t::ReinterpretCoords, g::Primitive) = TransformedGeometry(g, t)

applycoord(t::ReinterpretCoords, g::RegularGrid) = TransformedGrid(g, t)

applycoord(t::ReinterpretCoords, g::RectilinearGrid) = TransformedGrid(g, t)

applycoord(t::ReinterpretCoords, g::StructuredGrid) = TransformedGrid(g, t)

# -----------
# IO METHODS
# -----------

Base.show(io::IO, ::ReinterpretCoords{CRS₁,CRS₂}) where {CRS₁,CRS₂} =
  print(io, "ReinterpretCoords(CRS₁: $CRS₁, CRS₂: $CRS₂)")

function Base.show(io::IO, ::MIME"text/plain", t::ReinterpretCoords{CRS₁,CRS₂}) where {CRS₁,CRS₂}
  summary(io, t)
  println(io)
  println(io, "├─ CRS₁: $CRS₁")
  print(io, "└─ CRS₂: $CRS₂")
end
