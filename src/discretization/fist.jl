# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    FIST

Fast Industrial-Strength Triangulation (FIST) of polygons.

This triangulation method is the method behind the famous Mapbox's
Earcut library. It is based on a ear clipping algorithm adapted
for complex n-gons with holes. It has O(n²) time complexity where
n is the number of vertices. In practice it is very efficient due
to heuristics implemented in the algorithm.

## References

* Held. 1998. [FIST: Fast Industrial-Strength Triangulation of Polygons]
  (https://link.springer.com/article/10.1007/s00453-001-0028-4)
* Eder et al. 2018. [Parallelized ear clipping for the triangulation and
  constrained Delaunay triangulation of polygons]
  (https://www.sciencedirect.com/science/article/pii/S092577211830004X)
"""
struct FIST <: DiscretizationMethod end

function discretize(polygon::Polygon, method::FIST)
  verts, perms = _fist_remove_duplicates(polygon)
end

function _fist_remove_duplicates(polygon)
  outer, inners = rings(polygon)

  # retrieve vertices from rings
  verts = [@view vertices(r)[begin:end-1] for r in [outer; inners]]

  # sort vertices lexicographically
  perms = [sortperm(coordinates.(v)) for v in verts]

  # remove true duplicates
  vertsperms = map(1:length(verts)) do k
    vert = verts[k]
    perm = perms[k]

    keep = Int[]
    newperm = deepcopy(perm)
    sorted = @view vert[perm]
    for i in 1:length(sorted)-1
      if sorted[i] != sorted[i+1]
        # save index in the original vector
        push!(keep, perm[i])
      else
        # pop index from permutation vector
        # by assigning the value 0 and then
        # decrease all other indices that
        # are greater than the index by 1
        newperm[newperm .> perm[i]] .-= 1
        newperm[i] = 0
      end
    end
    push!(keep, last(perm))

    sort!(keep)
    filter!(!iszero, newperm)

    view(vert, keep), newperm
  end

  first.(vertsperms), last.(vertsperms)
end
