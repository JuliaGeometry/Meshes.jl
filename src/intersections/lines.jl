# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

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
    return @IT CrossingLines (a + λ₁ * (b - a)) f
  elseif r == rₐ == 1
    return @IT OverlappingLines line₁ f
  else
    return @IT NoIntersection nothing f
  end
end

# The intersection type can be one of six types:
# 1. intersect at one inner point (CrossingLineSegment -> Point)
# 2. intersect at an end point of segment (TouchingLineSegment -> Point)
# 3. overlap of line and segment (OverlappingLineSegment -> Segment)
# 4. do not overlap nor intersect (NoIntersection)
function intersection(f, line::Line{N,T}, seg::Segment{N,T}) where {N,T}
  a, b = line(0), line(1)
  c, d = seg(0), seg(1)

  # normalize points to gain parameter λ₂ corresponding to arc lengths
  l₂ = length(seg)
  d₀ = c + 1 / l₂ * (d - c)

  _, λ₂, r, rₐ = intersectparameters(a, b, c, d₀)

  # not in same plane or parallel
  if r ≠ rₐ
    return @IT NoIntersection nothing f # CASE 4
  # collinear
  elseif r == rₐ == 1
    return @IT OverlappingLineSegment seg f # CASE 3
  # in same plane, not parallel
  else
    λ₂ = mayberound(mayberound(λ₂, zero(T)), l₂)
    if λ₂ > 0 && λ₂ < l₂
      return @IT CrossingLineSegment seg(λ₂ / l₂) f # CASE 1, equal to line(λ₁)
    elseif λ₂ == 0 || λ₂ == l₂
      return @IT TouchingLineSegment ((λ₂ == 0) ? c : d) f # CASE 2
    else
      return @IT NoIntersection nothing f # CASE 4
    end
  end
end

# The intersection type can be one of four types:
# 
# 1. intersect at one inner point (CrossingRayLine -> Point)
# 2. intersect at origin of ray (TouchingRayLine -> Point)
# 3. overlap of line and ray (OverlappingRayLine -> Ray)
# 4. do not overlap nor intersect (NoIntersection)
function intersection(f, ray::Ray{N,T}, line::Line{N,T}) where {N,T}
  a, b = ray(0), ray(1)
  c, d = line(0), line(1)

  # rescaling of point b necessary to gain a parameter λ₁ representing the arc length
  l₁ = norm(b - a)
  b₀ = a + 1 / l₁ * (b - a)

  λ₁, _, r, rₐ = intersectparameters(a, b₀, c, d)

  if r ≠ rₐ # not in same plane or parallel
    return @IT NoIntersection nothing f # CASE 4
  elseif r == rₐ == 1 # collinear
    return @IT OverlappingRayLine ray f # CASE 3
  else # in same plane, not parallel
    λ₁ = mayberound(λ₁, zero(T))
    if λ₁ > 0
      return @IT CrossingRayLine ray(λ₁ / l₁) f # CASE 1
    elseif λ₁ == 0
      return @IT TouchingRayLine origin(ray) f # CASE 2
    else
      return @IT NoIntersection nothing f # CASE 4
    end
  end
end

# The intersection type can be one of five types:
#
# 1. intersect at one inner point (CrossingSegments -> Point)
# 2. intersect at one endpoint of one segment (MidTouchingSegments -> Point)
# 3. intersect at one endpoint of both segments (CornerTouchingSegments -> Point)
# 4. overlap of segments (OverlappingSegments -> Segments)
# 5. do not overlap nor intersect (NoIntersection)
function intersection(f, seg₁::Segment{N,T}, seg₂::Segment{N,T}) where {N,T}
  a, b = seg₁(0), seg₁(1)
  c, d = seg₂(0), seg₂(1)

  l₁, l₂ = length(seg₁), length(seg₂)
  b₀ = a + 1 / l₁ * (b - a) # corresponds to seg₁(1/length)
  d₀ = c + 1 / l₂ * (d - c)

  # arc length parameters λ₁ ∈ [0, l₁], λ₂ ∈ [0, l₂]: 
  λ₁, λ₂, r, rₐ = intersectparameters(a, b₀, c, d₀)

  if r ≠ rₐ # not in same plane or parallel
    return @IT NoIntersection nothing f #CASE 5
  elseif r == rₐ == 1 # collinear
    # find parameters λc and λd for points c and d in seg₁
    # use dimension with largest vector component to avoid division by zero
    v = b₀ - a
    i = argmax(abs.(v))
    λc, λd = ((c - a)[i], (d - a)[i]) ./ v[i]
    λc = mayberound(mayberound(λc, zero(T)), l₁)
    λd = mayberound(mayberound(λd, zero(T)), l₁)
    if (λc > l₁ && λd > l₁) || (λc < 0 && λd < 0)
      return @IT NoIntersection nothing f # CASE 5
    elseif (λc == 0 && λd < 0) || (λd == 0 && λc < 0)
      return @IT CornerTouchingSegments a f # CASE 3
    elseif (λc == l₁ && λd > l₁) || (λd == l₁ && λc > l₁)
      return @IT CornerTouchingSegments b f # CASE 3
    else
      params = sort([0, 1, λc / l₁, λd / l₁])
      p₁ = seg₁(params[2])
      p₂ = seg₁(params[3])
      return @IT OverlappingSegments Segment(p₁, p₂) f # CASE 4
    end
  else # in same plane, not parallel
    λ₁ = mayberound(mayberound(λ₁, zero(T)), l₁)
    λ₂ = mayberound(mayberound(λ₂, zero(T)), l₂)
    if λ₁ < 0 || λ₂ < 0 || λ₁ > l₁ || λ₂ > l₂
      return @IT NoIntersection nothing f # CASE 5
    # 8 cases remain
    elseif λ₁ == 0
      if λ₂ == 0 || λ₂ == l₂
        return @IT CornerTouchingSegments a f # CASE 3
      else
        return @IT MidTouchingSegments a f # CASE 2
      end
    elseif λ₁ == l₁
      if λ₂ == 0 || λ₂ == l₂
        return @IT CornerTouchingSegments b f # CASE 3
      else
        return @IT MidTouchingSegments b f # CASE 2
      end
    elseif λ₂ == 0 || λ₂ == l₂
      return @IT MidTouchingSegments (λ₂ == 0 ? c : d) f # CASE 2
    else
      return @IT CrossingSegments seg₁(λ₁ / l₁) f # CASE 1: equal to seg₂(λ₂/l₂)
    end
  end
end

# The intersection type can be one of six types:
#
# 1. intersect at one inner point (CrossingRays -> Point)
# 2. intersect at origin of one ray (MidTouchingRays -> Point)
# 3. intersect at origin of both rays (CornerTouchingRays -> Point)
# 4. overlap with aligned vectors (OverlappingAgreeingRays -> Ray)
# 5. overlap with colliding vectors (OverlappingOpposingRays -> Segment)
# 6. do not overlap nor intersect (NoIntersection)
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
    return @IT NoIntersection nothing f #CASE 6
  # collinear
  elseif r == rₐ == 1
    if direction(ray₁) ⋅ direction(ray₂) ≥ 0 # rays aligned in same direction
      if (origin(ray₁) - origin(ray₂)) ⋅ direction(ray₁) ≥ 0 # origin of ray₁ ∈ ray₂
        return @IT OverlappingAgreeingRays ray₁ f # CASE 4: ray₁
      else
        return @IT OverlappingAgreeingRays ray₂ f # CASE 4: ray₂
      end
    else # colliding rays
      if origin(ray₁) ∉ ray₂
        return @IT NoIntersection nothing f # CASE 6
      elseif origin(ray₁) == origin(ray₂)
        return @IT CornerTouchingRays a f # CASE 3
      else
        return @IT OverlappingOpposingRays Segment(origin(ray₁), origin(ray₂)) f # CASE 5
      end
    end
    # in same plane, not parallel
  else
    λ₁ = mayberound(λ₁, zero(T))
    λ₂ = mayberound(λ₂, zero(T))
    if λ₁ < 0 || λ₂ < 0
      return @IT NoIntersection nothing f # CASE 6
    elseif λ₁ == 0
      if λ₂ == 0
        return @IT CornerTouchingRays a f # CASE 3
      else
        return @IT MidTouchingRays a f # CASE 2: origin of ray₁
      end
    else
      if λ₂ == 0
        return @IT MidTouchingRays c f # CASE 2: origin of ray₂
      else
        return @IT CrossingRays ray₁(λ₁ / l₁) f # CASE 1: equal to ray₂(λ₂/l₂)
      end
    end
  end
end

# The intersection type can be one of five types:
# 
# 1. intersect at one inner point (CrossingRaySegment -> Point)
# 2. intersect at one corner point of segment xor origin of ray (MidTouchingRaySegment -> Point)
# 3. intersects at one corner point of segment and origin of ray (CornerTouchingRaySegment -> Point)
# 4. overlap at more than one point (OverlappingRaySegment -> Segment)
# 5. do not overlap nor intersect (NoIntersection)
function intersection(f, ray::Ray{N,T}, seg::Segment{N,T}) where {N,T}
  a, b = ray(0), ray(1)
  c, d = seg(0), seg(1)

  # normalize points to gain parameters λ₁, λ₂ corresponding to arc lengths
  l₁, l₂ = norm(b - a), length(seg)
  b₀ = a + 1 / l₁ * (b - a)
  d₀ = c + 1 / l₂ * (d - c)

  λ₁, λ₂, r, rₐ = intersectparameters(a, b₀, c, d₀)

  # not in same plane or parallel
  if r ≠ rₐ
    return @IT NoIntersection nothing f # CASE 5
  # collinear
  elseif r == rₐ == 1
    rc = sum((c - a) ./ direction(ray)) / N
    rd = sum((d - a) ./ direction(ray)) / N
    rc = mayberound(rc, zero(T))
    rd = mayberound(rd, zero(T))
    if rc > 0 # c ∈ ray
      if rd ≥ 0
        return @IT OverlappingRaySegment seg f # CASE 4
      else
        return @IT OverlappingRaySegment Segment(origin(ray), c) f # CASE 4
      end
    elseif rc == 0
      if rd > 0
        return @IT OverlappingRaySegment seg f # CASE 4
      else
        return @IT CornerTouchingRaySegment a f # CASE 3
      end
    else # rc < 0
      if rd > 0
        return @IT OverlappingRaySegment (Segment(origin(ray), d)) f # CASE 4
      elseif rd == 0
        return @IT CornerTouchingRaySegment a f # CASE 3
      else
        return @IT NoIntersection nothing f
      end
    end
    # in same plane, not parallel
  else
    λ₁ = mayberound(λ₁, zero(T))
    λ₂ = mayberound(mayberound(λ₂, zero(T)), l₂)
    if λ₁ < 0 || (λ₂ < 0 || λ₂ > l₂)
      return @IT NoIntersection nothing f
    elseif λ₁ == 0
      if λ₂ == 0 || λ₂ == l₂
        return @IT CornerTouchingRaySegment a f # CASE 3
      else
        return @IT MidTouchingRaySegment a f # CASE 2
      end
    else
      if λ₂ == 0 || λ₂ == l₂
        return @IT MidTouchingRaySegment (λ₂ < l₂ / 2 ? c : d) f # CASE 2
      else
        return @IT CrossingRaySegment ray(λ₁ / l₁) f # CASE 1, equal to seg(λ₂/l₂)
      end
    end
  end
end
