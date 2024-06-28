# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Planefy()

Convert the of coordinates of a geometry or domain to Cartesian2D when possible.
"""
struct Planefy <: CoordinateTransform end

applycoord(::Planefy, v::Vec) = v

applycoord(::Planefy, p::Point) = Point(_planefy(coords(p)))

_planefy(::CRS) = throw(ArgumentError("points must have 2 coordinates"))
_planefy(coords::Cartesian{Datum,2}) where {Datum} = coords
_planefy(coords::Polar) = convert(Cartesian, coords)
_planefy(coords::LatLon) = Cartesian{datum(coords)}(CoordRefSystems.rawvalues(coords))
_planefy(coords::GeocentricLatLon) = Cartesian{datum(coords)}(CoordRefSystems.rawvalues(coords))
_planefy(coords::AuthalicLatLon) = Cartesian{datum(coords)}(CoordRefSystems.rawvalues(coords))
_planefy(coords::CoordRefSystems.Projected) = convert(Cartesian, coords)
