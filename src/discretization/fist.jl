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

* Held, M. 1998. [FIST: Fast Industrial-Strength Triangulation of Polygons]
  (https://link.springer.com/article/10.1007/s00453-001-0028-4)
* Eder et al. 2018. [Parallelized ear clipping for the triangulation and
  constrained Delaunay triangulation of polygons]
  (https://www.sciencedirect.com/science/article/pii/S092577211830004X)
"""
struct FIST <: DiscretizationMethod end

function discretize(polyarea::PolyArea, ::FIST)
  # build bridges in case the polygonal area has
  # holes, i.e. reduce to a single outer boundary
  𝒫 = polyarea |> unique |> bridge

  # perform ear clipping
  while nvertices(𝒫) > 3
    # current vertices
    vs = vertices(𝒫)

    # CE1.1: classify angles as convex vs. reflex
    isconvex = innerangles(𝒫) .< π

    # CE1.2: check if segment vᵢ-₁ -- vᵢ+₁ intersects 𝒫
    intersects = map(1:nvertices(𝒫)) do i
      # target segment vᵢ-₁ -- vᵢ+₁
      sᵢ = Segment(vs[i-1], vs[i+1])

      # loop over all edges of 𝒫
      cross = false
      for j in 1:nvertices(𝒫)
        sⱼ = Segment(vs[j], vs[j+1])
        I = intersecttype(sᵢ, sⱼ)
        if !(I isa CornerTouchingSegments || I isa NonIntersectingSegments)
          cross = true
          break
        end
      end

      cross
    end
  end
end
