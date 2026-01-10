# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

include("polygonboolean/fillsegments.jl")
include("polygonboolean/utils.jl")
include("polygonboolean/sweep.jl")
include("polygonboolean/selection.jl")
include("polygonboolean/reconstruction.jl")

"""
    polygonbooleanop(geometry, other, operation)

Performs boolean `operation` of `geometry` with `other`.

## Operands

- `operation` — the operation to perform; one of:
  - `intersect` — regions inside both geometries (∩)
  - `union` — regions inside at least one geometry (∪)
  - `setdiff` — regions inside `geometry` but not `other`
  - `symdiff` — regions inside exactly one geometry (xor/⊻)

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

Base.union(a::Geometry, b::Geometry) = polygonbooleanop(a, b, union)
Base.setdiff(a::Geometry, b::Geometry) = polygonbooleanop(a, b, setdiff)
Base.symdiff(a::Geometry, b::Geometry) = polygonbooleanop(a, b, symdiff)
Base.xor(a::Geometry, b::Geometry) = polygonbooleanop(a, b, xor)

function polygonbooleanop(g1::Geometry, g2::Geometry, operation)
  g1 = deepcopy(g1)
  g2 = deepcopy(g2)
  polygonbooleanop!(g1, g2, operation)
end

function polygonbooleanop!(g1::Geometry, g2::Geometry, operation)
  throw(ArgumentError("$(operation) not supported for $(typeof(g1)) and $(typeof(g2))"))
end

# TODO: possible support for other being line segments that split poly
function polygonbooleanop!(poly::Polygon, other::Geometry, operation)
  subjrings = _getrings(poly)
  cliprings = _getrings(other)

  res = _buildrings(subjrings, cliprings, operation)
  isnothing(res) && return nothing

  segments, fills = res
  _buildpolys(segments, fills, operation)
end

function polygonbooleanop!(ring::Ring, other::Ring, operation)
  subjrings = _getrings(ring)
  cliprings = _getrings(other)

  res = _buildrings(subjrings, cliprings, operation)
  isnothing(res) && return nothing

  segments, fills = res
  polys = _buildpolys(segments, fills, operation)

  if polys isa PolyArea && length(rings(polys)) == 1
    return first(rings(polys))
  end

  polys
end

function _buildrings(subjrings, cliprings, operation)
  # flatten rings into segments for intersection
  subjsegs = [seg for r in subjrings for seg in Meshes.segments(Ring(r))]
  clipsegs = [seg for r in cliprings for seg in Meshes.segments(Ring(r))]

  # intersect
  allsegs = [subjsegs; clipsegs]
  intersections, seginds = pairwiseintersect(allsegs)

  # insert intersections into rings
  allrings = [subjrings; cliprings]
  _insertintersections!(intersections, seginds, allrings)

  # re-create segments from modified rings
  newsubjsegs = [seg for r in subjrings for seg in Meshes.segments(Ring(r))]
  newclipsegs = [seg for r in cliprings for seg in Meshes.segments(Ring(r))]

  nsegsfirst = length(newsubjsegs)
  newallsegs = [newsubjsegs; newclipsegs]

  # fill
  fillsegs = FillSegments(newallsegs)
  _annotatefill!(fillsegs, nsegsfirst)

  # select
  segments, fills = _selectsegments(fillsegs, operation)
  isempty(segments) && return nothing

  # remove degenerate segments
  keep = [s[1] != s[2] for s in segments]
  if !all(keep)
    segments = segments[keep]
    fills = fills[keep]
  end

  segments, fills
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
