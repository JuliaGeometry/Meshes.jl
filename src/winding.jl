# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    winding(points, object)

Generalized winding number of `points` with respect to the geometric `object`.

## References

* Barill et al. 2018. [Fast winding numbers for soups and clouds]
  (https://dl.acm.org/doi/10.1145/3197517.3201337)
"""
function winding end

function winding(p::Point{2,T}, r::Ring{2,T}) where {T}
  v = vertices(r)
  n = length(v)
  sum(∠(v[i], p, v[i + 1]) for i in 1:n) / T(2π)
end

# fallback with iterable of points
winding(points, object) = [winding(point, object) for point in points]
