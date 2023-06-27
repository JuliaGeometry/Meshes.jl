# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# The intersection type can be one of six types:
# 
# 1. Touching - origin of ray intersects middle of triangle
# 2. EdgeTouching - origin of ray intersects edge of triangle
# 3. CornerTouching - origin of ray intersects corner of triangle
# 4. Crossing - middle of ray intersects middle of triangle
# 5. EdgeCrossing - middle of ray intersects edge of triangle
# 6. CornerCrossing - middle of ray intersects corner of triangle
#
# Möller, T. & Trumbore, B., 1997.
# (https://www.tandfonline.com/doi/abs/10.1080/10867651.1997.10487468)
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
    return @IT NotIntersecting nothing f
  end

  # calculate u parameter and test bounds
  u = τ ⋅ p
  if u < -atol(T) || u > det
    return @IT NotIntersecting nothing f
  end

  q = τ × e₁

  # calculate v parameter and test bounds
  v = d ⋅ q
  if v < -atol(T) || u + v > det
    return @IT NotIntersecting nothing f
  end

  λ = (e₂ ⋅ q) * (one(T) / det)

  if λ < -atol(T)
    return @IT NotIntersecting nothing f
  end

  # assemble barycentric weights
  w = Vec(u, v, det - u - v)

  if any(isapprox.(o, vs, atol=atol(T)))
    return @IT CornerTouching r(λ) f
  elseif isapprox(λ, zero(T), atol=atol(T))
    if all(>(zero(T)), w)
      return @IT Touching r(λ) f
    else
      return @IT EdgeTouching r(λ) f
    end
  end

  if count(x -> isapprox(x, zero(T), atol=atol(T)), w) == 1
    return @IT EdgeCrossing r(λ) f
  elseif count(x -> isapprox(x, det, atol=atol(T)), w) == 1
    return @IT CornerCrossing r(λ) f
  end

  λ = clamp(λ, zero(T), typemax(T))

  return @IT Crossing r(λ) f
end

# Jiménez, J., Segura, R. and Feito, F. 2009.
# (https://www.sciencedirect.com/science/article/pii/S0925772109001448?via%3Dihub)
function intersection(f, s::Segment{3,T}, t::Triangle{3,T}) where {T}
  vₛ = vertices(s)
  vₜ = vertices(t)

  A = vₛ[1] - vₜ[3]
  B = vₜ[1] - vₜ[3]
  C = vₜ[2] - vₜ[3]

  W₁ = B × C
  wᵥ = A ⋅ W₁

  D = vₛ[2] - vₜ[3]
  sᵥ = D ⋅ W₁

  if wᵥ > atol(T)
    # rejection 2
    if sᵥ > atol(T)
      return @IT NotIntersecting nothing f
    end

    W₂ = A × D
    tᵥ = W₂ ⋅ C

    # rejection 3
    if tᵥ < -atol(T)
      return @IT NotIntersecting nothing f
    end

    uᵥ = -(W₂ ⋅ B)

    # rejection 4
    if uᵥ < -atol(T)
      return @IT NotIntersecting nothing f
    end

    # rejection 5
    if wᵥ < (sᵥ + tᵥ + uᵥ)
      return @IT NotIntersecting nothing f
    end
  elseif wᵥ < -atol(T)
    # rejection 2
    if sᵥ < -atol(T)
      return @IT NotIntersecting nothing f
    end

    W₂ = A × D
    tᵥ = W₂ ⋅ C

    # rejection 3
    if tᵥ > atol(T)
      return @IT NotIntersecting nothing f
    end

    uᵥ = -(W₂ ⋅ B)

    # rejection 4
    if uᵥ > atol(T)
      return @IT NotIntersecting nothing f
    end

    # rejection 5
    if wᵥ > (sᵥ + tᵥ + uᵥ)
      return @IT NotIntersecting nothing f
    end
  else
    if sᵥ > atol(T)
      W₂ = D × A
      tᵥ = W₂ ⋅ C

      # rejection 3
      if tᵥ < -atol(T)
        return @IT NotIntersecting nothing f
      end

      uᵥ = -(W₂ ⋅ B)

      # rejection 4
      if uᵥ < -atol(T)
        return @IT NotIntersecting nothing f
      end
      # rejection 5
      if -sᵥ < (tᵥ + uᵥ)
        return @IT NotIntersecting nothing f
      end
    elseif sᵥ < -atol(T)
      W₂ = D × A
      tᵥ = W₂ ⋅ C

      # rejection 3
      if tᵥ > atol(T)
        return @IT NotIntersecting nothing f
      end

      uᵥ = -(W₂ ⋅ B)

      # rejection 4
      if uᵥ > atol(T)
        return @IT NotIntersecting nothing f
      end

      # rejection 5
      if -sᵥ > (tᵥ + uᵥ)
        return @IT NotIntersecting nothing f
      end
    else
      # rejection 1, coplanar segment
      return @IT NotIntersecting nothing f
    end
  end

  λ = clamp(wᵥ / (wᵥ - sᵥ), zero(T), one(T))

  return @IT Intersecting s(λ) f
end
