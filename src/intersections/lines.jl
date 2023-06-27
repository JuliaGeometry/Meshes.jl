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
  a, b = seg₁(0), seg₁(1)
  c, d = seg₂(0), seg₂(1)

  l₁, l₂ = length(seg₁), length(seg₂)
  b₀ = a + 1 / l₁ * (b - a) # corresponds to seg₁(1/length)
  d₀ = c + 1 / l₂ * (d - c)

  # arc length parameters λ₁ ∈ [0, l₁], λ₂ ∈ [0, l₂]: 
  λ₁, λ₂, r, rₐ = intersectparameters(a, b₀, c, d₀)

  if r ≠ rₐ # not in same plane or parallel
    return @IT NotIntersecting nothing f #CASE 5
  elseif r == rₐ == 1 # collinear
    # find parameters λc and λd for points c and d in seg₁
    # use dimension with largest vector component to avoid division by zero
    v = b₀ - a
    i = argmax(abs.(v))
    λc, λd = ((c - a)[i], (d - a)[i]) ./ v[i]
    λc = mayberound(mayberound(λc, zero(T)), l₁)
    λd = mayberound(mayberound(λd, zero(T)), l₁)
    if (λc > l₁ && λd > l₁) || (λc < 0 && λd < 0)
      return @IT NotIntersecting nothing f # CASE 5
    elseif (λc == 0 && λd < 0) || (λd == 0 && λc < 0)
      return @IT CornerTouching a f # CASE 3
    elseif (λc == l₁ && λd > l₁) || (λd == l₁ && λc > l₁)
      return @IT CornerTouching b f # CASE 3
    else
      params = sort([0, 1, λc / l₁, λd / l₁])
      p₁ = seg₁(params[2])
      p₂ = seg₁(params[3])
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
#
# 1. intersect at one inner point (Crossing -> Point)
# 2. intersect at origin of one ray (EdgeTouching -> Point)
# 3. intersect at origin of both rays (CornerTouching -> Point)
# 4. overlap with aligned vectors (PosOverlapping -> Ray)
# 5. overlap with colliding vectors (NegOverlapping -> Segment)
# 6. do not overlap nor intersect (NotIntersecting -> Nothing)
function intersection(f, ray₁::Ray{N,T}, ray₂::Ray{N,T}) where {N,T}
  a, b = ray₁(0), ray₁(1)
  c, d = ray₂(0), ray₂(1)

  # normalize points to gain parameters λ₁, λ₂ corresponding to arc lengths
  l₁, l₂ = norm(b - a), norm(d - c)
  b₀ = a + 1 / l₁ * (b - a)
  d₀ = c + 1 / l₂ * (d - c)

  λ₁, λ₂, r, rₐ = intersectparameters(a, b₀, c, d₀)

  # not in same plane or parallel
  if r ≠ rₐ
    return @IT NotIntersecting nothing f #CASE 6
  # collinear
  elseif r == rₐ == 1
    if (b - a) ⋅ (d - c) ≥ 0 # rays aligned in same direction
      if (a - c) ⋅ (b - a) ≥ 0 # origin of ray₁ ∈ ray₂
        return @IT PosOverlapping ray₁ f # CASE 4: ray₁
      else
        return @IT PosOverlapping ray₂ f # CASE 4: ray₂
      end
    else # colliding rays
      if a ∉ ray₂
        return @IT NotIntersecting nothing f # CASE 6
      elseif a == c
        return @IT CornerTouching a f # CASE 3
      else
        return @IT NegOverlapping Segment(a, c) f # CASE 5
      end
    end
    # in same plane, not parallel
  else
    λ₁ = mayberound(λ₁, zero(T))
    λ₂ = mayberound(λ₂, zero(T))
    if λ₁ < 0 || λ₂ < 0
      return @IT NotIntersecting nothing f # CASE 6
    elseif λ₁ == 0
      if λ₂ == 0
        return @IT CornerTouching a f # CASE 3
      else
        return @IT EdgeTouching a f # CASE 2: origin of ray₁
      end
    else
      if λ₂ == 0
        return @IT EdgeTouching c f # CASE 2: origin of ray₂
      else
        return @IT Crossing ray₁(λ₁ / l₁) f # CASE 1: equal to ray₂(λ₂/l₂)
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

# The intersection type can be one of four types:
# 
# 1. intersect at one inner point (Crossing -> Point)
# 2. intersect at origin of ray (Touching -> Point)
# 3. overlap of line and ray (Overlapping -> Ray)
# 4. do not overlap nor intersect (NotIntersecting -> Nothing)
function intersection(f, ray::Ray{N,T}, line::Line{N,T}) where {N,T}
  a, b = ray(0), ray(1)
  c, d = line(0), line(1)

  # rescaling of point b necessary to gain a parameter λ₁ representing the arc length
  l₁ = norm(b - a)
  b₀ = a + 1 / l₁ * (b - a)

  λ₁, _, r, rₐ = intersectparameters(a, b₀, c, d)

  if r ≠ rₐ # not in same plane or parallel
    return @IT NotIntersecting nothing f # CASE 4
  elseif r == rₐ == 1 # collinear
    return @IT Overlapping ray f # CASE 3
  else # in same plane, not parallel
    λ₁ = mayberound(λ₁, zero(T))
    if λ₁ > 0
      return @IT Crossing ray(λ₁ / l₁) f # CASE 1
    elseif λ₁ == 0
      return @IT Touching a f # CASE 2
    else
      return @IT NotIntersecting nothing f # CASE 4
    end
  end
end

# The intersection type can be one of three types:
# 
# 1. intersect at one point
# 2. overlap at more than one point
# 3. do not overlap nor intersect
function intersection(f, line₁::Line, line₂::Line)
  a, b = line₁(0), line₁(1)
  c, d = line₂(0), line₂(1)

  λ₁, _, r, rₐ = intersectparameters(a, b, c, d)

  if r == rₐ == 2
    return @IT Crossing (a + λ₁ * (b - a)) f
  elseif r == rₐ == 1
    return @IT Overlapping line₁ f
  else
    return @IT NotIntersecting nothing f
  end
end
