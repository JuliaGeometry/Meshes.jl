# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    withcrs(g, coords, CRS=Cartesian)

Point with the same CRS of `g` from another point with `coords` in given `CRS`.
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
  _fromvalues(first(points), values)
end

"""
    coordmean(points; weights=nothing)
  
Mean of the base coordinates of the points, `Cartesian` for `𝔼` and `LatLon` for `🌐`.
If `weights` is passed, the weighted mean will be returned.
"""
function coordmean(points; weights=nothing)
  den = if isnothing(weights)
    length(points)
  else
    sum(weights)
  end
  values = _coordsum(points, weights) ./ den
  _fromvalues(first(points), values)
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

function _tovalues(p)
  CRS = _basecrs(manifold(p))
  c = convert(CRS, coords(p))
  CoordRefSystems.values(c)
end

function _fromvalues(g, values)
  CRS = _basecrs(manifold(g))
  withcrs(g, values, CRS)
end

function _coordsum(points, weights)
  if isnothing(weights)
    mapreduce(_tovalues, .+, points)
  else
    mapreduce((p, w) -> _tovalues(p) .* w, .+, points, weights)
  end
end

_basecrs(::Type{<:𝔼}) = Cartesian
_basecrs(::Type{<:🌐}) = LatLon
