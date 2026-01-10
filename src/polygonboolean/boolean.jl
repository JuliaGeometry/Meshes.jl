# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
  PolygonBoolean

Types for geometric boolean operations.

  * PolyIntersection - regions inside both geometries
  * PolyUnion - regions inside at least one geometry
  * PolyDifference - regions inside geometry but not other
  * PolySymDifference - regions inside exactly one geometry
"""
abstract type PolygonBoolean end
struct PolyIntersection <: PolygonBoolean end
struct PolyUnion <: PolygonBoolean end
struct PolyDifference <: PolygonBoolean end
struct PolySymDifference <: PolygonBoolean end

_retrieveoperation(::PolyIntersection) = :intersection
_retrieveoperation(::PolyUnion) = :union
_retrieveoperation(::PolyDifference) = :difference
_retrieveoperation(::PolySymDifference) = :xor

"""
    polygonbooleanop(geometry, other, operation)

Performs boolean `operation` from `geometry` with `other`.

## Operands

- `operation` — the operation to perform; one of:
  - `PolyIntersection` — regions inside both geometries (intersect/∩)
  - `PolyUnion` — regions inside at least one geometry (union/∪)
  - `PolyDifference` — regions inside `geometry` but not `other` (set difference/setdiff)
  - `PolySymDifference` — regions inside exactly one geometry (symmetric difference/xor/⊻)

## Notes

The algorithm works for both convex and concave polygons using
Martinez-Rueda clipping. For mutating version, use `polygonbooleanop!`.

## References
* Martínez, F., Rueda, A.J., Feito, F.R. 2009. [A new algorithm for computing Boolean operations on
  polygons](https://doi.org/10.1016/j.cag.2009.03.003)
"""
function polygonbooleanop end

"""
  g₁ ∪ g₂

Compute the geometric union of geometries `g₁` and `g₂`.
"""
Base.union

"""
  setdiff(g₁, g₂)

Compute the geometric difference of geometries `g₁` and `g₂`.
"""
Base.setdiff

"""
  symdiff(g₁, g₂)

Compute the symmetric difference of geometries `g₁` and `g₂`.
"""
Base.symdiff

"""
  g₁ ⊻ g₂

Compute the symmetric difference (exclusive or) of geometries `g₁` and `g₂`.
"""
Base.xor

Base.union(a::Geometry, b::Geometry) = polygonbooleanop(a, b, PolyUnion())
Base.setdiff(a::Geometry, b::Geometry) = polygonbooleanop(a, b, PolyDifference())
Base.symdiff(a::Geometry, b::Geometry) = polygonbooleanop(a, b, PolySymDifference())
Base.xor(a::Geometry, b::Geometry) = polygonbooleanop(a, b, PolySymDifference())

function polygonbooleanop(g1::Geometry, g2::Geometry, operation::PolygonBoolean)
  g1 = deepcopy(g1)
  g2 = deepcopy(g2)
  polygonbooleanop!(g1, g2, operation)
end

function polygonbooleanop!(g1::Geometry, g2::Geometry, operation::PolygonBoolean)
  throw(ArgumentError("$(operation) not supported for $(typeof(g1)) and $(typeof(g2))"))
end

function polygonbooleanop!(poly::Polygon, other::Geometry, operation::PolygonBoolean)
  # Extract rings
  subjrings = _getrings(poly)
  cliprings = _getrings(other)

  # Flatten rings into segments for intersection
  subjsegs = [Segment(r[i], r[i < length(r) ? i + 1 : 1]) for r in subjrings for i in 1:length(r)]
  clipsegs = [Segment(r[i], r[i < length(r) ? i + 1 : 1]) for r in cliprings for i in 1:length(r)]

  # Intersect
  allsegs = [subjsegs; clipsegs]
  intersections, seginds = pairwiseintersect(allsegs)

  # Insert intersections into rings
  allrings = [subjrings; cliprings]
  _insertintersections!(intersections, seginds, allrings)

  # Re-create segments from modified rings
  newsubjsegs = [Segment(r[i], r[i < length(r) ? i + 1 : 1]) for r in subjrings for i in 1:length(r)]
  newclipsegs = [Segment(r[i], r[i < length(r) ? i + 1 : 1]) for r in cliprings for i in 1:length(r)]

  nsubj = length(newsubjsegs)
  newallsegs = [newsubjsegs; newclipsegs]

  # Fill
  fillsegs = FillSegments(newallsegs)
  _annotatefill!(fillsegs, nsubj)

  # Select
  segments, fills, _ = _selectsegments(fillsegs, nsubj, operation)
  isempty(segments) && return nothing

  # Remove degenerate segments
  keep = [s[1] != s[2] for s in segments]
  if !all(keep)
    segments = segments[keep]
    fills = fills[keep]
  end

  _buildrings(segments, fills, operation)
end

function polygonbooleanop!(ring::Ring, other::Ring, operation::PolygonBoolean)
  polygonbooleanop!(PolyArea(ring), PolyArea(other), operation)
end

function _getrings(poly::Polygon)
  [collect(vertices(r)) for r in rings(poly)]
end

function _getrings(geom::Geometry)
  _getrings(convert(Quadrangle, geom))
end

function _getrings(ring::Ring)
  [collect(vertices(ring))]
end
