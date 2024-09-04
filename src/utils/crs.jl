# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    withcrs(g, coords, basecrs=Cartesian)

Point with the same CRS of `g` from another point with `coords` specified in a source `basecrs`.
"""
function withcrs(g::GeometryOrDomain, coords::Tuple; basecrs=Cartesian)
  M = manifold(g)
  C = crs(g)
  D = datum(C)
  c = convert(C, basecrs{D}(coords...))
  Point{M}(c)
end

"""
    withcrs(g, v)

Point at the end of the vector `v` with the same CRS of `g`.
"""
withcrs(g::GeometryOrDomain, v::StaticVector) = withcrs(g, Tuple(v), basecrs=Cartesian)

"""
    flat(p)

Flatten coordinates of point `p` to Cartesian coordinates,
ignoring the original units of the coordinate reference system.
"""
flat(p::Point) = Point(flat(coords(p)))
flat(c::CRS) = Cartesian{datum(c)}(CoordRefSystems.raw(c))
