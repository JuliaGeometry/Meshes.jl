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

  # points of resulting mesh
  points = vertices(𝒫)

  # perform ear clipping
  𝒬 = ears(𝒫)
  𝒯 = Connectivity{Triangle,3}[]
  while nvertices(𝒫) > 3
    clipped = false
    if !isempty(𝒬)
      i = pop!(𝒬)
      push!(𝒯, connect((i-1,i,i+1), Triangle))
      clipped = true
    elseif clipped
      𝒬 = ears(𝒫)
    else
      # recovery process
      @warn "entered in recovery process"
    end
  end
end

function ears(𝒫)
  v = vertices(𝒫)

  # CE1.1: classify angles as convex vs. reflex
  isconvex = innerangles(𝒫) .< π

  # CE1.2: check if segment vᵢ-₁ -- vᵢ+₁ intersects 𝒫
  intersects = map(1:nvertices(𝒫)) do i
    sᵢ = Segment(v[i-1], v[i+1])
    cross = false
    for j in 1:nvertices(𝒫)
      sⱼ = Segment(v[j], v[j+1])
      I = intersecttype(sᵢ, sⱼ)
      if !(I isa CornerTouchingSegments || I isa NonIntersectingSegments)
        cross = true
        break
      end
    end
    cross
  end

  # CE1.3: check if vᵢ-1 ∈ C(vᵢ, vᵢ+1, vᵢ+2) and vᵢ+1 ∈ C(vᵢ-2, vᵢ-1, vᵢ)
  incone = map(1:nvertices(𝒫)) do i
    c1 = sideof(v[i-1], Segment(v[i+1], v[i  ])) != :LEFT
    c2 = sideof(v[i-1], Segment(v[i+1], v[i+2])) != :RIGHT
    c3 = sideof(v[i+1], Segment(v[i-1], v[i-2])) != :LEFT
    c4 = sideof(v[i+1], Segment(v[i-1], v[i  ])) != :RIGHT
    all((c1, c2, c3, c4))
  end

  findall(isconvex .& .!intersects .& incone)
end
