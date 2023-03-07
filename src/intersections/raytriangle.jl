# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

#=
Möller, T. & Trumbore, B., 1997.
(https://www.tandfonline.com/doi/abs/10.1080/10867651.1997.10487468)

Cases
1. CrossingRayTriangle - middle of ray intersects middle of triangle
2. TouchingRayTriangle - origin of ray intersects middle of triangle
3. CornerTouchingRayTriangle - origin of ray intersects corner of triangle
4. EdgeTouchingRayTriangle - origin of ray intersects edge of triangle
5. EdgeCrossingRayTriangle - middle of ray intersects edge of triangle
6. CornerCrossingRayTriangle - middle of ray intersects corner of triangle
=#
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

  # assemble barycentric weights
  w = Vec(u, v, det - u - v)

  if any(isapprox.(o, vs, atol=atol(T)))
    return @IT CornerTouchingRayTriangle r(λ) f
  elseif isapprox(λ, zero(T), atol=atol(T))
    if all(w .> zero(T))
      return @IT TouchingRayTriangle r(λ) f
    else
      return @IT EdgeTouchingRayTriangle r(λ) f
    end
  end

  if count(x -> isapprox(x, zero(T), atol=atol(T)), w) == 1
    return @IT EdgeCrossingRayTriangle r(λ) f
  elseif count(x -> isapprox(x, det, atol=atol(T)), w) == 1
    return @IT CornerCrossingRayTriangle r(λ) f
  end

  λ = clamp(λ, zero(T), typemax(T))

  return @IT CrossingRayTriangle r(λ) f
end
