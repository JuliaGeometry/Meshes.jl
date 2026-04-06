# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    JarvisMarch()
    JarvisMarch(k)

Compute the convex hull of a set of points or geometries using the
Jarvis's march algorithm. See [https://en.wikipedia.org/wiki/Gift_wrapping_algorithm]
(https://en.wikipedia.org/wiki/Gift_wrapping_algorithm).

If `k` is provided, the algorithm will attempt to compute a concave hull using [k,n)
nearest neighbors.

The algorithm has complexity `O(n*h)` where `n` is the number of points
and `h` is the number of points in the hull.

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

struct AdaptiveJarvisMarch{K<:JarvisMarch{<:Integer}} <: HullMethod
  march::K
end
AdaptiveJarvisMarch(k::I) where {I<:Integer} = AdaptiveJarvisMarch(JarvisMarch(k))

function hull(points, method::AdaptiveJarvisMarch)
  k = method.march.k
  kmax = length(points)
  while true
    try
      chul = hull(points, JarvisMarch(k))
      all(p .∈ Ref(chul) for p in points) && return chul
      throw(ArgumentError("Not all points are in hull with k = $k. Try increasing k."))
    catch e
      if e isa ArgumentError
        k += 1
        k > kmax && return hull(points, JarvisMarch())
      else
        rethrow(e)
      end
    end
  end
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
  i₀ = copy(i)
  # initialize hull with i
  ℐ = [i]

  # initialize candidate set and searcher for k-nearest neighbors if needed
  k = method.k
  if k isa Integer
    T = typeof(k)
    k = min(max(k, 3), n)
    assertion(ispositive(k), "k must be a positive integer")
    k > T(n - 2) && return hull(points, JarvisMarch()) # fallback to convex hull
    searcher = KNearestSearch(p, k)
    𝒞 = trues(n)
  else
    searcher = nothing
    𝒞 = 1:n
  end
  candidates = jarviscandidates(searcher, 𝒞, n, p, i)
  # find next point with smallest angle
  O = p[i]
  A = O + Vec(zero(ℒ), -oneunit(ℒ))
  j = jarvisnext(searcher, candidates, p, ℐ, A, O, k)

  # initialize ring of indices
  push!(ℐ, j)

  # rotational sweep
  step = 0 # only used for k-nearest approach but low cost inclusion
  while first(ℐ) != last(ℐ)
    # if k-nearest neighbors are being used, reinsert first point after a bit
    k isa Integer && step == 5 && (𝒞[i₀] = true)
    # direction of current segment
    v = p[j] - p[i]

    # update candidates
    candidates = jarviscandidates(searcher, 𝒞, n, p, j)
    isempty(candidates) && throw(
      ArgumentError(
        "Not enough points to compute hull with k = $k. Try increasing k or computing the convex hull instead."
      )
    )

    # find next segment
    i = j
    O = p[i]
    A = O + v
    j = jarvisnext(searcher, candidates, p, ℐ, A, O, k)

    # update ring of indices
    push!(ℐ, j)
    step += 1
  end

  # return polygonal area
  PolyArea(p[ℐ[begin:(end - 1)]])
end

# helper to find next point for concave hull
function jarvisnext(::KNearestSearch, candidates, p, ℐ, A, O, k)
  sort!(candidates, by=l -> ∠(A, O, p[l]))
  valid = nothing
  # if fewer than 3 points, we can skip intersection checks
  length(ℐ) < 3 && return first(candidates)
  # check candidates in order of angle until we find one that doesn't intersect the existing hull
  for idx in eachindex(candidates)
    nᵢ = candidates[idx]
    cpoint = p[nᵢ]
    offset = cpoint == p[ℐ[begin]] ? 1 : 0
    limit = length(ℐ) - 1 - offset
    ok = limit < 2 || !any(2:limit) do indⱼ
      intersects(Segment(p[ℐ[end]], cpoint), Segment(p[ℐ[end - indⱼ + 1]], p[ℐ[end - indⱼ]]))
    end
    if ok
      valid = idx
      break
    end
  end
  isnothing(valid) && throw(ArgumentError("No valid candidate point found for hull construction with k = $k"))
  candidates[valid]
end

jarvisnext(::Nothing, candidates, p, ℐ, A, O, k) = argmin(l -> ∠(A, O, p[l]), candidates)

# helper to get candidate indices for next point
jarviscandidates(searcher::Nothing, 𝒞, n, p, i) = setdiff(1:n, i)

function jarviscandidates(searcher::KNearestSearch, 𝒞, n, p, i)
  # mask current point to avoid returning it as a candidate
  𝒞[i] = false # updates 𝒞 in place so prior points lead to early termination
  search(p[i], searcher; mask=𝒞)
end
