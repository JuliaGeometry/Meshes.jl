# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FlatCoords()

Convert the of coordinates of a geometry or domain to Cartesian2D when possible.
"""
struct FlatCoords <: CoordinateTransform end

applycoord(::FlatCoords, v::Vec) = v

applycoord(::FlatCoords, p::Point) = Point(_flatcoords(coords(p)))

_flatcoords(::CRS) = throw(ArgumentError("points must have 2 coordinates"))
_flatcoords(coords::Cartesian{Datum,2}) where {Datum} = coords
_flatcoords(coords::Polar) = convert(Cartesian, coords)
_flatcoords(coords::LatLon) = Cartesian{datum(coords)}(CoordRefSystems.rawvalues(coords))
_flatcoords(coords::GeocentricLatLon) = Cartesian{datum(coords)}(CoordRefSystems.rawvalues(coords))
_flatcoords(coords::AuthalicLatLon) = Cartesian{datum(coords)}(CoordRefSystems.rawvalues(coords))
_flatcoords(coords::CoordRefSystems.Projected) = convert(Cartesian, coords)
