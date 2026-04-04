# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    withcrs(g, coords, CRS=basecrs(g))

Point with the same CRS of `g` from another point with `coords` in given `CRS`.
"""
withcrs(g::GeometryOrDomain, coords::Tuple, CRS=basecrs(g)) = Point{manifold(g)}(convert(crs(g), CRS(coords...)))

"""
    withcrs(g, v)

Point at the end of the vector `v` with the same CRS of `g`.
"""
withcrs(g::GeometryOrDomain, v::StaticVector) = withcrs(g, Tuple(v), Cartesian{datum(crs(g))})

"""
    basecrs(g)

Base coordinate reference system of `g` as a function of the manifold.
"""
function basecrs(p::Point)
  D = datum(crs(p))
  manifold(p) === 🌐 ? LatLon{D} : Cartesian{D}
end

basecrs(g::Geometry) = basecrs(centroid(g))

basecrs(d::Domain) = basecrs(first(d))

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
  coordvals(p) = CoordRefSystems.values(convert(basecrs(p), coords(p)))
  if isnothing(weights)
    mapreduce(coordvals, .+, points)
  else
    mapreduce((p, w) -> coordvals(p) .* w, .+, points, weights)
  end
end
