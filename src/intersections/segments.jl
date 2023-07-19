# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# The intersection type can be one of five types:
#
# 1. intersect at one inner point (Crossing -> Point)
# 2. intersect at one endpoint of one segment (EdgeTouching -> Point)
# 3. intersect at one endpoint of both segments (CornerTouching -> Point)
# 4. overlap of segments (Overlapping -> Segments)
# 5. do not overlap nor intersect (NotIntersecting -> Nothing)
function intersection(f, seg₁::Segment{N,T}, seg₂::Segment{N,T}) where {N,T}
  a, b = vertices(seg₁)
  c, d = vertices(seg₂)

  # handle degenerate segments
  if a == b && c == d
    if a ≈ c
      return @IT CornerTouching a f
    else
      return @IT NotIntersecting nothing f
    end
  elseif a == b
    if a ≈ c || a ≈ d
      return @IT CornerTouching a f
    elseif a ∈ seg₂
      return @IT EdgeTouching a f
    else
      return @IT NotIntersecting nothing f
    end
  elseif c == d
    if c ≈ a || c ≈ b
      return @IT CornerTouching c f
    elseif c ∈ seg₁
      return @IT EdgeTouching c f
    else
      return @IT NotIntersecting nothing f
    end
  end

  l₁ = length(seg₁)
  l₂ = length(seg₂)
  b₀ = a + 1 / l₁ * (b - a)
  d₀ = c + 1 / l₂ * (d - c)

  # arc length parameters λ₁ ∈ [0, l₁], λ₂ ∈ [0, l₂]: 
  λ₁, λ₂, r, rₐ = intersectparameters(a, b₀, c, d₀)

  if r ≠ rₐ # not in same plane or parallel
    return @IT NotIntersecting nothing f #CASE 5
  elseif r == rₐ == 1 # collinear
    # find parameters λc and λd for points c and d in seg₁
    # use dimension with largest vector component to avoid division by zero
    v = b₀ - a
    i = last(findmax(abs, v))
    vc = c - a
    vd = d - a
    λc = vc[i] / v[i]
    λd = vd[i] / v[i]
    λc = mayberound(mayberound(λc, zero(T)), l₁)
    λd = mayberound(mayberound(λd, zero(T)), l₁)
    if (λc > l₁ && λd > l₁) || (λc < 0 && λd < 0)
      return @IT NotIntersecting nothing f # CASE 5
    elseif (λc == 0 && λd < 0) || (λd == 0 && λc < 0)
      return @IT CornerTouching a f # CASE 3
    elseif (λc == l₁ && λd > l₁) || (λd == l₁ && λc > l₁)
      return @IT CornerTouching b f # CASE 3
    else
      t₁, t₂ = _sort4vals(zero(T), one(T), λc / l₁, λd / l₁)
      p₁ = seg₁(t₁)
      p₂ = seg₁(t₂)
      return @IT Overlapping Segment(p₁, p₂) f # CASE 4
    end
  else # in same plane, not parallel
    λ₁ = mayberound(mayberound(λ₁, zero(T)), l₁)
    λ₂ = mayberound(mayberound(λ₂, zero(T)), l₂)
    if λ₁ < 0 || λ₂ < 0 || λ₁ > l₁ || λ₂ > l₂
      return @IT NotIntersecting nothing f # CASE 5
    # 8 cases remain
    elseif λ₁ == 0
      if λ₂ == 0 || λ₂ == l₂
        return @IT CornerTouching a f # CASE 3
      else
        return @IT EdgeTouching a f # CASE 2
      end
    elseif λ₁ == l₁
      if λ₂ == 0 || λ₂ == l₂
        return @IT CornerTouching b f # CASE 3
      else
        return @IT EdgeTouching b f # CASE 2
      end
    elseif λ₂ == 0 || λ₂ == l₂
      return @IT EdgeTouching (λ₂ == 0 ? c : d) f # CASE 2
    else
      return @IT Crossing seg₁(λ₁ / l₁) f # CASE 1: equal to seg₂(λ₂/l₂)
    end
  end
end

# The intersection type can be one of five types:
# 
# 1. intersect at one inner point (Crossing -> Point)
# 2. intersect at one end point of segment xor origin of ray (EdgeTouching -> Point)
# 3. intersects at one end point of segment and origin of ray (CornerTouching -> Point)
# 4. overlap at more than one point (Overlapping -> Segment)
# 5. do not overlap nor intersect (NotIntersecting -> Nothing)
function intersection(f, seg::Segment{N,T}, ray::Ray{N,T}) where {N,T}
  a, b = ray(0), ray(1)
  c, d = seg(0), seg(1)

  # normalize points to gain parameters λ₁, λ₂ corresponding to arc lengths
  l₁, l₂ = norm(b - a), length(seg)
  b₀ = a + 1 / l₁ * (b - a)
  d₀ = c + 1 / l₂ * (d - c)

  λ₁, λ₂, r, rₐ = intersectparameters(a, b₀, c, d₀)

  # not in same plane or parallel
  if r ≠ rₐ
    return @IT NotIntersecting nothing f # CASE 5
  # collinear
  elseif r == rₐ == 1
    rc = sum((c - a) ./ (b - a)) / N
    rd = sum((d - a) ./ (b - a)) / N
    rc = mayberound(rc, zero(T))
    rd = mayberound(rd, zero(T))
    if rc > 0 # c ∈ ray
      if rd ≥ 0
        return @IT Overlapping seg f # CASE 4
      else
        return @IT Overlapping Segment(ray(0), c) f # CASE 4
      end
    elseif rc == 0
      if rd > 0
        return @IT Overlapping seg f # CASE 4
      else
        return @IT CornerTouching a f # CASE 3
      end
    else # rc < 0
      if rd > 0
        return @IT Overlapping (Segment(ray(0), d)) f # CASE 4
      elseif rd == 0
        return @IT CornerTouching a f # CASE 3
      else
        return @IT NotIntersecting nothing f
      end
    end
    # in same plane, not parallel
  else
    λ₁ = mayberound(λ₁, zero(T))
    λ₂ = mayberound(mayberound(λ₂, zero(T)), l₂)
    if λ₁ < 0 || (λ₂ < 0 || λ₂ > l₂)
      return @IT NotIntersecting nothing f
    elseif λ₁ == 0
      if λ₂ == 0 || λ₂ == l₂
        return @IT CornerTouching a f # CASE 3
      else
        return @IT EdgeTouching a f # CASE 2
      end
    else
      if λ₂ == 0 || λ₂ == l₂
        return @IT EdgeTouching (λ₂ < l₂ / 2 ? c : d) f # CASE 2
      else
        return @IT Crossing ray(λ₁ / l₁) f # CASE 1, equal to seg(λ₂/l₂)
      end
    end
  end
end

# The intersection type can be one of six types:
# 1. intersect at one inner point (Crossing -> Point)
# 2. intersect at an end point of segment (Touching -> Point)
# 3. overlap of line and segment (Overlapping -> Segment)
# 4. do not overlap nor intersect (NotIntersecting -> Nothing)
function intersection(f, seg::Segment{N,T}, line::Line{N,T}) where {N,T}
  a, b = line(0), line(1)
  c, d = seg(0), seg(1)

  # normalize points to gain parameter λ₂ corresponding to arc lengths
  l₂ = length(seg)
  d₀ = c + 1 / l₂ * (d - c)

  _, λ₂, r, rₐ = intersectparameters(a, b, c, d₀)

  # not in same plane or parallel
  if r ≠ rₐ
    return @IT NotIntersecting nothing f # CASE 4
  # collinear
  elseif r == rₐ == 1
    return @IT Overlapping seg f # CASE 3
  # in same plane, not parallel
  else
    λ₂ = mayberound(mayberound(λ₂, zero(T)), l₂)
    if λ₂ > 0 && λ₂ < l₂
      return @IT Crossing seg(λ₂ / l₂) f # CASE 1, equal to line(λ₁)
    elseif λ₂ == 0 || λ₂ == l₂
      return @IT Touching ((λ₂ == 0) ? c : d) f # CASE 2
    else
      return @IT NotIntersecting nothing f # CASE 4
    end
  end
end

# Jiménez, J., Segura, R. and Feito, F. 2009.
# (https://www.sciencedirect.com/science/article/pii/S0925772109001448?via%3Dihub)
function intersection(f, seg::Segment{3,T}, tri::Triangle{3,T}) where {T}
  vₛ = vertices(seg)
  vₜ = vertices(tri)

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

  return @IT Intersecting seg(λ) f
end

# sorts four numbers using a sorting network 
# and returns the 2nd and 3rd
function _sort4vals(a, b, c, d)
  a > c && ((a, c) = (c, a))
  b > d && ((b, d) = (d, b))
  a > b && ((a, b) = (b, a))
  c > d && ((c, d) = (d, c))
  b > c && ((b, c) = (c, b))
  b, c
end
