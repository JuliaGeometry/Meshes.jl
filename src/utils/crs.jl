# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    withcrs(g, coords, CRS=Cartesian)

Point with the same CRS of `g` from another point with `coords` specified in `CRS`.
"""
function withcrs(g::GeometryOrDomain, coords::Tuple, ::Type{CRS}) where {CRS}
  M = manifold(g)
  C = crs(g)
  D = datum(C)
  c = convert(C, CRS{D}(coords...))
  Point{M}(c)
end

withcrs(g::GeometryOrDomain, coords::Tuple) = withcrs(g, coords, Cartesian)

"""
    withcrs(g, v)

Point at the end of the vector `v` with the same CRS of `g`.
"""
withcrs(g::GeometryOrDomain, v::StaticVector) = withcrs(g, Tuple(v), Cartesian)

"""
    flat(p)

Flatten coordinates of point `p` to Cartesian coordinates,
ignoring the original units of the coordinate reference system.
"""
flat(p::Point) = Point(flat(coords(p)))
flat(c::CRS) = Cartesian{datum(c)}(CoordRefSystems.raw(c))
