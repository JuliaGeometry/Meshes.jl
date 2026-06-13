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
nearest neighbors. However, the algorithm is not guaranteed to succeed for any particular `k`.

The algorithm has complexity `O(n*h)` where `n` is the number of points
and `h` is the number of points in the hull.

see `concavehull` for a version that iteratively increases `k` until all
points are in the hull, which is useful when an effective `k` is not known prior.
This gurantees correctness.

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
  kₒ = method.k
  isnothing(kₒ) || kₒ ≥ n && return hull(points, JarvisMarch())

  # find bottom-left point
  i = argmin(p)
  start = i
  # initialize hull with i
  ℐ = [i]

  # initialize mask of candidate points and searcher for k-nearest neighbors if needed
  k, searcher, pointmask = jarvissearcher(method.k, n, p)

  # find neighbor candidates
  𝒞 = jarviscandidates(searcher, pointmask, n, p, i, i, start)

  # find next point with smallest angle
  O = p[i]
  A = O + Vec(zero(ℒ), -oneunit(ℒ))
  j = jarvisnext(searcher, 𝒞, p, ℐ, A, O, k)

  # initialize ring of indices
  push!(ℐ, j)
  k isa Integer && (pointmask[j] = true)

  # rotational sweep
  while first(ℐ) != last(ℐ)
    # direction of current segment
    v = p[j] - p[i]

    # find candidates for next point
    inds = isnothing(searcher) ? (i, j) : j
    𝒞 = jarviscandidates(searcher, pointmask, n, p, j, inds, start)
    isempty(𝒞) && return nothing # no candidates, should only happen if k is too small

    # find next segment
    i = j
    O = p[i]
    A = O + v
    j = jarvisnext(searcher, 𝒞, p, ℐ, A, O, k)
    isnothing(j) && return nothing # no valid next point, should only happen if k is too small

    # update ring of indices
    push!(ℐ, j)
    k isa Integer && (pointmask[j] = true)
  end

  # return polygonal area
  PolyArea(p[ℐ[begin:(end - 1)]])
end

# helper to find next point for convex hull
jarvisnext(::Nothing, 𝒞, p, ℐ, A, O, k) = argmin(l -> ∠(A, O, p[l]), 𝒞)

function jarvisnext(::KNearestSearch, 𝒞, p, ℐ, A, O, k)
  sort!(𝒞, by=l -> ∠(A, O, p[l]))
  # check candidates in order of angle until we find one that doesn't intersect the existing hull
  for nᵢ in 𝒞
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
jarviscandidates(::Nothing, pointmask, n, p, current, inds, start) = setdiff(1:n, inds)

function jarviscandidates(searcher::KNearestSearch, pointmask, n, p, current, inds, start)
  # mask out points already in the hull, except for the starting point
  mask = .!pointmask
  mask[start] = true
  mask[inds] = false
  search(p[current], searcher; mask=mask)
end

# outputs k, searcher, and a mask of candidate point indices
jarvissearcher(k, n, p) = nothing, nothing, 1:n

function jarvissearcher(kₒ::T, n, p) where {T<:Integer}
  k = min(max(kₒ, T(3)), n)
  assertion(ispositive(k), "k must be a positive integer")
  k, KNearestSearch(p, k), falses(n)
end
