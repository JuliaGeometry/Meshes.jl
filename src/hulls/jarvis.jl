# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    JarvisMarch()
    JarvisMarch(k)

Compute the convex hull of a set of points or geometries using the
Jarvis's march algorithm. See [https://en.wikipedia.org/wiki/Gift_wrapping_algorithm]
(https://en.wikipedia.org/wiki/Gift_wrapping_algorithm).

If `k` is provided, the algorithm will attempt to compute a concave hull using k
nearest neighbors.

The algorithm has complexity `O(n*h)` where `n` is the number of points
and `h` is the number of points in the hull.

see `AdaptiveJarvisMarch` for a version that iteratively increases `k` until all
points are in the hull, which is useful when an effective `k` is not known prior.

## References

* Jarvis 1973. [On the identification of the convex hull of a finite set of
  points in the plane](https://www.sciencedirect.com/science/article/abs/pii/0020019073900203)
* Moreira, A. & Santos, M. Y. 2007. "concave hull: a k-nearest neighbours approach for the computation of the region occupied by a set of points". In: *Proceedings of the Second International Conference on Computer Graphics Theory and Applications*. SciTePress, pp. 61–68.
"""
struct JarvisMarch{K} <: HullMethod
  k::K
end
JarvisMarch() = JarvisMarch{Nothing}(nothing)
JarvisMarch(k::I) where {I<:Integer} = JarvisMarch{I}(k)

"""
    AdaptiveJarvisMarch(k)

Compute the concave hull of a set of points or geometries using an adaptive
version of the Jarvis's march algorithm. The algorithm will iteratively increase `k` until all points are in the hull

see also `JarvisMarch` for creating a `convexhull` and `concavehull` for default
  concave hulls starting with `k = 3`.

"""
struct AdaptiveJarvisMarch{K<:JarvisMarch{<:Integer}} <: HullMethod
  march::K
end
AdaptiveJarvisMarch() = AdaptiveJarvisMarch(JarvisMarch(3))
AdaptiveJarvisMarch(k::I) where {I<:Integer} = AdaptiveJarvisMarch(JarvisMarch(k))

function hull(points, method::AdaptiveJarvisMarch)
  k = method.march.k
  pointsꜝ = unique(points)
  for kᵢ in k:length(pointsꜝ)
    chul = hull(points, JarvisMarch(kᵢ))
    # validate we found a hull that contains all points, if so return it, otherwise increase k and try again
    isnothing(chul) && continue
    all(points .∈ Ref(chul)) && return chul
  end
  # otherwise, fallback to convex hull
  hull(points, JarvisMarch())
end

function hull(points, method::JarvisMarch)
  pₒ = first(points)
  ℒ = lentype(pₒ)

  # sanity check
  ncoords = CoordRefSystems.ncoords(coords(pₒ))
  assertion(ncoords == 2, "Jarvis's march algorithm is only defined with 2D coordinates")

  # remove duplicates
  p = unique(points)
  n = length(p)

  # corner cases
  n == 1 && return p[1]
  n == 2 && return Segment(p[1], p[2])
  n == 3 && return PolyArea(p)

  # find bottom-left point
  i = argmin(p)
  start = i
  # initialize hull with i
  ℐ = [i]

  # initialize candidate set and searcher for k-nearest neighbors if needed
  k, searcher, 𝒞 = _jarvissetup(method.k, n, p)
  # if no candidates exist (k too large), fallback to convex hull
  isnothing(𝒞) && return hull(points, JarvisMarch())

  candidates = jarviscandidates(searcher, 𝒞, n, p, i, start)

  # find next point with smallest angle
  O = p[i]
  A = O + Vec(zero(ℒ), -oneunit(ℒ))
  j = jarvisnext(searcher, candidates, p, ℐ, A, O, k)

  # initialize ring of indices
  push!(ℐ, j)
  k isa Integer && (𝒞[j] = true)

  # rotational sweep
  while first(ℐ) != last(ℐ)
    # direction of current segment
    v = p[j] - p[i]

    # find candidates for next point
    candidates = jarviscandidates(searcher, 𝒞, n, p, j, start)
    isempty(candidates) && return nothing # no candidates, should only happen if k is too small

    # find next segment
    i = j
    O = p[i]
    A = O + v
    j = jarvisnext(searcher, candidates, p, ℐ, A, O, k)
    isnothing(j) && return nothing # no valid next point, should only happen if k is too small

    # update ring of indices
    push!(ℐ, j)
    k isa Integer && (𝒞[j] = true)
  end

  # return polygonal area
  PolyArea(p[ℐ[begin:(end - 1)]])
end

# helper to find next point for concave hull
jarvisnext(::Nothing, candidates, p, ℐ, A, O, k) = argmin(l -> ∠(A, O, p[l]), candidates)

function jarvisnext(::KNearestSearch, candidates, p, ℐ, A, O, k)
  sort!(candidates, by=l -> ∠(A, O, p[l]))
  # check candidates in order of angle until we find one that doesn't intersect the existing hull
  for nᵢ in candidates
    cpoint = p[nᵢ]
    cseg = Segment(p[ℐ[end]], cpoint)
    cbox = boundingbox(cseg)
    offset = cpoint == p[ℐ[begin]] ? 1 : 0
    limit = length(ℐ) - 1 - offset
    ok = limit < 2 || !any(2:limit) do indⱼ
      p₁ = p[ℐ[end - indⱼ + 1]]
      p₂ = p[ℐ[end - indⱼ]]
      eseg = Segment(p₁, p₂)
      # quick check to see if segments could intersect before doing more expensive segment intersection check
      intersects(cbox, boundingbox(eseg)) || return false
      intersects(cseg, eseg)
    end
    ok && return nᵢ
  end
  nothing
end

# helper to get candidate indices for next point
jarviscandidates(::Nothing, 𝒞, n, p, i, start) = setdiff(1:n, i)

function jarviscandidates(searcher::KNearestSearch, 𝒞, n, p, i, start)
  # mask out points already in the hull, except for the starting point
  mask = .!𝒞
  mask[start] = true
  mask[i] = false
  search(p[i], searcher; mask=mask)
end

# outputs k, searcher, and a mask of candidate point indices
_jarvissetup(k, n, p) = nothing, nothing, 1:n

function _jarvissetup(k::T, n, p) where {T<:Integer}
  k = min(max(k, 3), n)
  assertion(ispositive(k), "k must be a positive integer")
  k, KNearestSearch(p, k), falses(n)
end
