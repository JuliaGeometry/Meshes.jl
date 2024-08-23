# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    withcrs(g, v)

Point at the end of the vector `v` with the same CRS of `g`.
"""
function withcrs(g::GeometryOrDomain, v::StaticVector)
  C = crs(g)
  D = datum(C)
  c = convert(C, Cartesian{D}(v...))
  Point(c)
end

"""
    flat(p)

Flatten coordinates of point `p` to Cartesian coordinates,
ignoring the original units of the coordinate reference system.
"""
flat(p::Point) = Point(flat(coords(p)))
flat(c::CRS) = Cartesian{datum(c)}(CoordRefSystems.raw(c))
