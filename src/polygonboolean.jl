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
Base.xor(a::Geometry, b::Geometry) = polygonbooleanop(a, b, symdiff)

function polygonbooleanop(g1::Geometry, g2::Geometry, operation)
  g1 = deepcopy(g1)
  g2 = deepcopy(g2)
  polygonbooleanop!(g1, g2, operation)
end

function polygonbooleanop!(g1::Geometry, g2::Geometry, operation)
  throw(ArgumentError("$(operation) not supported for $(typeof(g1)) and $(typeof(g2))"))
end

# TODO: build up functions progressively over PRs
# TODO: possible support for other being line segments that split poly
# TODO: polygon vs ring method dispatch needs to be tweaked. best for now, both operands same type
function polygonbooleanop!(poly::Polygon, other::Geometry, operation)
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
  subjsegs = Meshes.segments.(subjrings)
  clipsegs = Meshes.segments.(cliprings)

  # intersect
  allsegs = Iterators.flatten((subjsegs, clipsegs))
  intersections, seginds = pairwiseintersect(allsegs)

  # insert intersections into rings
  allrings = vertices.((subjrings, cliprings))
  _insertintersections!(intersections, seginds, allrings)

  # annotate segments with fill information

  # cleanup of annotations
end
