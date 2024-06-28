# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FlatCoords()

Convert the coordinates of a geometry or domain to Cartesian2D when possible.
"""
struct FlatCoords <: CoordinateTransform end

applycoord(::FlatCoords, v::Vec) = v

applycoord(::FlatCoords, p::Point) = Point(_flatcoords(coords(p)))

function _flatcoords(coords::CRS) 
  if CoordRefSystems.ncoords(coords) ≠ 2
    throw(ArgumentError("points must have 2 coordinates"))
  end
  convert(Cartesian, coords)
end

_flatcoords(coords::LatLon) = Cartesian{datum(coords)}(CoordRefSystems.rawvalues(coords))

_flatcoords(coords::GeocentricLatLon) = Cartesian{datum(coords)}(CoordRefSystems.rawvalues(coords))

_flatcoords(coords::AuthalicLatLon) = Cartesian{datum(coords)}(CoordRefSystems.rawvalues(coords))
