# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    IntersectionType

The different types of sides that a point may lie in relation to a
boundary geometry or mesh. Type `SideType` in a Julia session to see
the full list.
"""
@enum SideType begin
  IN
  OUT
  ON
  LEFT
  RIGHT
end

"""
    sideof(points, object)

Determines on which side the `points` are in relation to the geometric
`object`, which can be a boundary `geometry` or `mesh`.
"""
function sideof end

# ---------
# GEOMETRY
# ---------

"""
    sideof(point, line)

Determines on which side the `point` is in relation to the `line`.
Possible results are `LEFT`, `RIGHT` or `ON` the `line`.

### Notes

* Assumes the orientation of `Segment(line(0), line(1))`.
"""
function sideof(point::Point{2,T}, line::Line{2,T}) where {T}
  a = signarea(point, line(0), line(1))
  ifelse(a > atol(T), LEFT, ifelse(a < -atol(T), RIGHT, ON))
end

"""
    sideof(point, ring)

Determines on which side the `point` is in relation to the `ring`.
Possible results are `IN` or `OUT` the `ring`.
"""
function sideof(point::Point{2,T}, ring::Ring{2,T}) where {T}
  w = winding(point, ring)
  ifelse(isapprox(w, zero(T), atol=atol(T)), OUT, IN)
end

"""
    sideof(point, simplex)

Determines on which side the `point` is in relation to a k-simplex.
Possible results are `LEFT` or `RIGHT`, although they are only useful for
comparison, not for geometrical interpretation.
"""
function sideof(point::Point{Dim}, splx::Simplex{K,Dim}) where {K,Dim}
  Dim == (K + 1) || throw(ArgumentError("`sideof` is only defined for `paramdim(splx)+1==Dim`"))
  ifelse(signbit(normal(splx)' * (centroid(splx) - point)), LEFT, RIGHT)
end

# -----
# MESH
# -----

"""
    sideof(point, mesh)

Determines on which side the `point` is in relation to the surface `mesh`.
Possible results are `IN` or `OUT` the `mesh`.
"""
sideof(point::Point{3}, mesh::Mesh{3}) = sideof((point,), mesh) |> first

# ----------
# FALLBACKS
# ----------

sideof(points, line::Line{2}) = map(point -> sideof(point, line), points)

function sideof(points, object::GeometryOrDomain)
  T = coordtype(object)
  bbox = boundingbox(object)
  isin = tcollect(point ∈ bbox for point in points)
  inds = findall(isin)
  wind = winding(collectat(points, inds), object)
  side = fill(OUT, length(isin))
  side[inds] .= ifelse.(isapprox.(wind, zero(T), atol=atol(T)), OUT, IN)
  side
end
