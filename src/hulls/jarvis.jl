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
  # find next point with smallest angle
  O = p[i]
  A = O + Vec(zero(ℒ), -oneunit(ℒ))

  # perform primary Jarvis march loop
  ℐ = _jarvisloop(method, points, p, n, i, O, A)
  # if no hull is found, return convex hull as fallback
  isnothing(ℐ) && return convexhull(p)

  # return polygonal area
  PolyArea(p[ℐ[begin:(end - 1)]])
end

# convex hull
function _jarvisloop(::JarvisMarch{Nothing}, points, p, n, i, O, A)
  # candidates for next point
  𝒞 = [1:(i - 1); (i + 1):n]

  # find next point with smallest angle
  j = argmin(l -> ∠(A, O, p[l]), 𝒞)

  # initialize ring of indices
  ℐ = [i, j]

  # rotational sweep
  while first(ℐ) != last(ℐ)
    # direction of current segment
    v = p[j] - p[i]

    # update candidates
    𝒞 = setdiff(1:n, [i, j])

    # find next segment
    i = j
    O = p[i]
    A = O + v
    j = argmin(l -> ∠(A, O, p[l]), 𝒞)

    # update ring of indices
    push!(ℐ, j)
  end

  # return indices of hull vertices
  ℐ
end

# concave hull
function _jarvisloop(method::JarvisMarch{I}, points, p, n, i, O, A) where {I<:Integer}
  m = I(n - 2)
  k = min(max(method.k, 3), m)
  assertion(ispositive(k), "k must be a positive integer")
  k > m && return _jarvisloop(JarvisMarch{Nothing}(nothing), points, p, n, i, O, A) # fallback to convex hull

  # initial state for retries
  i₀, O₀, A₀ = i, O, A

  # try increasing k until valid hull found
  for ki in k:m
    searcher = KNearestSearch(p, ki)
    mask = trues(n)

    # apply initial point information
    i = i₀
    O = O₀
    A = A₀
    mask[i] = false

    # find next point with smallest angle among k nearest neighbors
    𝒩 = search(O, searcher; mask=mask)
    isempty(𝒩) && (k += 1; continue)

    j = argmin(l -> ∠(A, O, p[l]), 𝒩)
    ℐ = [i, j]

    # no longer searchable
    mask[j] = false
    step = 2
    failed = false

    while first(ℐ) != last(ℐ)
      # reinsert first point after 5 steps. Otherwise the searcher may double back and fail to find a valid hull.
      step == 5 && (mask[ℐ[begin]] = true)

      # direction of current segment
      v = p[j] - p[i]
      i = j
      O = p[i]
      A = O + v

      # order neighbors by angle
      search!(𝒩, O, searcher; mask)
      isempty(𝒩) && break
      sort!(𝒩, by=l -> ∠(A, O, p[l]))

      found = false

      for nᵢ in 𝒩
        cpoint = p[nᵢ]
        last = cpoint == p[ℐ[begin]] ? 1 : 0

        its = false
        indⱼ = 2
        while !its && indⱼ < length(ℐ) - last
          its = intersects(Segment(p[ℐ[step]], cpoint), Segment(p[ℐ[step - indⱼ + 1]], p[ℐ[step - indⱼ]]))
          indⱼ += 1
        end

        if !its
          j = nᵢ
          found = true
          break
        end
      end

      # if no candidate found, break and try again with larger k
      !found && (failed = true; break)

      mask[j] = false
      step += 1
      push!(ℐ, j)
    end

    # check if hull is valid
    if !failed
      poly = PolyArea(p[ℐ[begin:(end - 1)]])
      all(points .∈ poly) && return ℐ
    end
  end

  nothing
end
