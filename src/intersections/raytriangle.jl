# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    intersection(f, ray, triangle)

Compute the intersection of a `ray` and `triangle`
and apply function `f` to it.

## References

* Möller, T. and Trumbore, B., 1997. [Fast, minimum storage ray-triangle
  intersection. Journal of graphics tools]
  (https://www.tandfonline.com/doi/abs/10.1080/10867651.1997.10487468)
"""
function intersection(f, r::Ray{3,T}, t::Triangle{3,T}) where {T}
  vs = vertices(t)
  o = origin(r)
  d = direction(r)

  e₁ = vs[3] - vs[1]
  e₂ = vs[2] - vs[1]
  p = d × e₂
  det = e₁ ⋅ p

  # keep det > 0, modify T accordingly
  if det > atol(T)
    τ = o - vs[1]
  else
    τ = vs[1] - o
    det = -det
  end

  if det < atol(T)
    # This ray is parallel to the plane of the triangle.
    return @IT NoIntersection nothing f
  end

  # calculate u parameter and test bounds
  u = τ ⋅ p
  if u < -atol(T) || u > det
    return @IT NoIntersection nothing f
  end

  q = τ × e₁

  # calculate v parameter and test bounds
  v = d ⋅ q
  if v < -atol(T) || u + v > det
    return @IT NoIntersection nothing f
  end

  λ = (e₂ ⋅ q) * (one(T) / det)

  if λ < -atol(T)
    return @IT NoIntersection nothing f
  end

  λ = clamp(λ, zero(T), typemax(T))

  return @IT IntersectingRayTriangle r(λ) f
end
