# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    withcrs(g, c, srccrs=manifoldcrs(g))

Point with the same CRS of geometry `g` converted from a
tuple of coordinates `c` in a given source CRS `srccrs`.
"""
function withcrs(g::GeometryOrDomain, c::Tuple, srccrs=manifoldcrs(g))
  dstcrs = basecrs(g)
  coords = convert(dstcrs, srccrs(c...))
  Point{manifold(g)}(coords)
end

"""
    withcrs(g, v)

Point at the end of the vector `v` with the same CRS of `g`.
"""
withcrs(g::GeometryOrDomain, v::StaticVector) = withcrs(g, Tuple(v), Cartesian{datum(crs(g))})

"""
    manifoldcrs(g)

Coordinate reference system for the manifold of geometry `g`.
"""
function manifoldcrs(g::GeometryOrDomain)
  D = datum(crs(g))
  manifold(g) === 🌐 ? LatLon{D} : Cartesian{D}
end

"""
    basecrs(g)

Coordinate reference system of geometry `g` without number type.
"""
basecrs(g::GeometryOrDomain) = CoordRefSystems.constructor(crs(g))

"""
    flat(p)

Flatten coordinates of point `p` to Cartesian coordinates,
ignoring the original units of the coordinate reference system.
"""
flat(p::Point) = Point(flat(coords(p)))
flat(c::CRS) = Cartesian{datum(c)}(CoordRefSystems.raw(c))

"""
    svec(p)

Return `SVector` with raw coordinates of point `p`.

### Notes

This utility function exists because NearestNeighbors.jl
currently only accepts coordinates of type `AbstractVector`.
"""
svec(p::Point) = svec(coords(p))
svec(c::CRS) = SVector(CoordRefSystems.raw(c))

"""
    coordsum(points; weights=nothing)
  
Sum of the base coordinates of the points, `Cartesian` for `𝔼` and `LatLon` for `🌐`.
If `weights` is passed, the weighted sum will be returned.
"""
function coordsum(points; weights=nothing)
  values = _coordsum(points, weights)
  withcrs(first(points), values)
end

"""
    coordmean(points; weights=nothing)
  
Mean of the base coordinates of the points, `Cartesian` for `𝔼` and `LatLon` for `🌐`.
If `weights` is passed, the weighted mean will be returned.
"""
function coordmean(points; weights=nothing)
  denom = isnothing(weights) ? length(points) : sum(weights)
  values = _coordsum(points, weights) ./ denom
  withcrs(first(points), values)
end

"""
    coordround(point, r=RoundNearest; digits=0, base=10)
    coordround(point, r=RoundNearest; sigdigits=0)

Round the coordinates of a `point` to specified presicion.
"""
function coordround(point::Point, r::RoundingMode=RoundNearest; kwargs...)
  c = coords(point)
  x = CoordRefSystems.values(c)
  x′ = round.(eltype(x), x, r; kwargs...)
  c′ = CoordRefSystems.constructor(c)(x′...)
  Point(c′)
end

# -----------------
# HELPER FUNCTIONS
# -----------------

function _coordsum(points, weights)
  coordvals(p) = CoordRefSystems.values(convert(manifoldcrs(p), coords(p)))
  if isnothing(weights)
    mapreduce(coordvals, .+, points)
  else
    mapreduce((p, w) -> coordvals(p) .* w, .+, points, weights)
  end
end
