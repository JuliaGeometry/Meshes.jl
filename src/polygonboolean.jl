# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# TODO: include subfolders for modularity
include("polygonboolean/utils.jl")
# include("polygonboolean/fillsegments.jl")
# include("polygonboolean/sweep.jl")
# include("polygonboolean/selection.jl")
# include("polygonboolean/reconstruction.jl")

"""
    polygonbooleanop(poly, other, operation)

Performs polygon boolean `operation` of `poly` with `other`.

## Operands

- `operation` — the operation to perform; one of:
  - `intersect` — regions inside both geometries (∩)
  - `union` — regions inside at least one geometry (∪)
  - `setdiff` — regions inside `poly` but not `other`
  - `symdiff` — regions inside exactly one geometry (xor/⊻)

## Notes

The algorithm works for both convex and concave polygons using
Martinez-Rueda clipping. For mutating version, use `polygonbooleanop!`.

## References
* Martínez, F., Rueda, A.J., Feito, F.R. 2009. [A new algorithm for computing Boolean operations on
  polygons](https://doi.org/10.1016/j.cag.2009.03.003)
"""
function polygonbooleanop end

function polygonbooleanop(g1::Polygon, g2::Polygon, operation)
  g1 = deepcopy(g1)
  g2 = deepcopy(g2)
  polygonbooleanop!(g1, g2, operation)
end

# fallback
function polygonbooleanop!(g1::Geometry, g2::Geometry, operation)
  throw(ArgumentError("$(operation) not supported for $(typeof(g1)) and $(typeof(g2))"))
end

# TODO: build up functions progressively over PRs
# TODO: possible support for `other` being line segments that split poly
# TODO: polygon vs ring method dispatch needs to be tweaked. best for now, both operands same type
function polygonbooleanop!(poly::Polygon, other::Polygon, operation)
  subjrings = rings(poly)
  cliprings = rings(other)

  res = _buildrings(subjrings, cliprings, operation)
end

function polygonbooleanop!(ring::Ring, other::Ring, operation)
  subjrings = rings(ring)
  cliprings = rings(other)

  # build rings with polygon relationship information
  res = _buildrings(subjrings, cliprings, operation)

  # build polygons from rings
end

function _buildrings(subjrings, cliprings, operation)
  # flatten rings into segments for intersection
  allrings = [subjrings; cliprings]
  allsegs = Iterators.flatten(segments.(allrings))
  intersections, seginds = pairwiseintersect(allsegs)

  # insert intersections into rings (mutates allrings)
  _insertintersections!(allrings, intersections, seginds)

  # annotate segments with fill information

  # cleanup of annotations
end
