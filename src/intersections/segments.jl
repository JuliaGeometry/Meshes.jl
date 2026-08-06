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
function intersection(f, seg₁::Segment, seg₂::Segment)
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

  l₁ = ustrip(length(seg₁))
  l₂ = ustrip(length(seg₂))
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
    λc = mayberound(mayberound(λc, zero(λc)), l₁)
    λd = mayberound(mayberound(λd, zero(λd)), l₁)
    if (λc > l₁ && λd > l₁) || (λc < 0 && λd < 0)
      return @IT NotIntersecting nothing f # CASE 5
    elseif (λc == 0 && λd < 0) || (λd == 0 && λc < 0)
      return @IT CornerTouching a f # CASE 3
    elseif (λc == l₁ && λd > l₁) || (λd == l₁ && λc > l₁)
      return @IT CornerTouching b f # CASE 3
    else
      t₁, t₂ = _sort4vals(zero(λc), one(λc), λc / l₁, λd / l₁)
      p₁ = seg₁(t₁)
      p₂ = seg₁(t₂)
      return @IT Overlapping Segment(p₁, p₂) f # CASE 4
    end
  else # in same plane, not parallel
    λ₁ = mayberound(mayberound(λ₁, zero(λ₁)), l₁)
    λ₂ = mayberound(mayberound(λ₂, zero(λ₂)), l₂)
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
function intersection(f, seg::Segment, ray::Ray)
  Dim = embeddim(seg)
  a, b = ray(0), ray(1)
  c, d = seg(0), seg(1)

  # normalize points to gain parameters λ₁, λ₂ corresponding to arc lengths
  l₁ = ustrip(norm(b - a))
  l₂ = ustrip(length(seg))
  b₀ = a + 1 / l₁ * (b - a)
  d₀ = c + 1 / l₂ * (d - c)

  λ₁, λ₂, r, rₐ = intersectparameters(a, b₀, c, d₀)

  # not in same plane or parallel
  if r ≠ rₐ
    return @IT NotIntersecting nothing f # CASE 5
  # collinear
  elseif r == rₐ == 1
    rc = sum((c - a) ./ (b - a)) / Dim
    rd = sum((d - a) ./ (b - a)) / Dim
    rc = mayberound(rc, zero(rc))
    rd = mayberound(rd, zero(rd))
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
    λ₁ = mayberound(λ₁, zero(λ₁))
    λ₂ = mayberound(mayberound(λ₂, zero(λ₂)), l₂)
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
function intersection(f, seg::Segment, line::Line)
  a, b = line(0), line(1)
  c, d = seg(0), seg(1)

  # normalize points to gain parameter λ₂ corresponding to arc lengths
  l₂ = ustrip(length(seg))
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
    λ₂ = mayberound(mayberound(λ₂, zero(λ₂)), l₂)
    if λ₂ > 0 && λ₂ < l₂
      return @IT Crossing seg(λ₂ / l₂) f # CASE 1, equal to line(λ₁)
    elseif λ₂ == 0 || λ₂ == l₂
      return @IT Touching ((λ₂ == 0) ? c : d) f # CASE 2
    else
      return @IT NotIntersecting nothing f # CASE 4
    end
  end
end

# intersection between a segment and a polygon
# 1. overlap with the polygon interior (Intersecting -> Segment or Multi)
# 2. overlap with a polygon edge only (EdgeTouching -> Segment or Multi)
# 3. intersect at an interior point of the segment and a polygon vertex (CornerTouching -> Point)
# 4. intersect at an endpoint of the segment (Touching -> Point)
# 5. the segment is degenerate and is inside the polygon (Touching -> Point)
# 6. do not overlap nor intersect (NotIntersecting -> Nothing)
function intersection(f, seg::Segment{𝔼{2}}, poly::Polygon{𝔼{2}})
  a, b = vertices(seg)
  segbbox = boundingbox(seg)
  polybbox = boundingbox(poly)

  # check for degenerate segments
  if a ≈ b && a ∈ poly
    return @IT Touching a f # Case 5
  end

  # assess obvious no intersection case beforehand
  if !intersects(segbbox, polybbox)
    return @IT NotIntersecting nothing f # Case 6
  end

  # segment endpoints as scalars + segment vector
  segλs = [0.0, 1.0]
  v = b - a
  v² = v ⋅ v

  # extract the first Cartesian coordinate
  x(p) = flat(coords(p)).x
  y(p) = flat(coords(p)).y

  # choose the polygon's longest bounding-box axis
  polymin, polymax = extrema(polybbox)
  xside = x(polymax) - x(polymin)
  yside = y(polymax) - y(polymin)
  axiscoord = xside ≥ yside ? x : y

  # segment range along the chosen axis
  segmin, segmax = extrema(segbbox)
  segstart = axiscoord(segmin)
  segstop = axiscoord(segmax)

  # assess intersections
  boundarypoints = typeof(a)[]
  boundaryoverlaps = Tuple{Float64,Float64}[]

  for ring in rings(poly)
    ringbbox = boundingbox(ring)

    # check for obvious no intersection case beforehand
    intersects(segbbox, ringbbox) || continue

    edges = collect(segments(ring))
    boxes = boundingbox.(edges)

    # edge ranges
    starts = map(boxes) do box
      pmin, _ = extrema(box)
      axiscoord(pmin)
    end
    stops = map(boxes) do box
      _, pmax = extrema(box)
      axiscoord(pmax)
    end

    # Sort by minimum coordinate
    inds = sortperm(starts)
    edges = edges[inds]
    boxes = boxes[inds]
    starts = starts[inds]
    stops = stops[inds]

    for i in eachindex(edges)
      # edge ends before the segment begins. Later edges may still overlap, so skip only this edge.
      stops[i] < segstart && continue

      # edge starts after the segment ends. Since edges are sorted by the selected axis minimum, all later edges are also too far right.
      starts[i] > segstop && break

      # the selected axis ranges overlap, but the full bounding boxes may still be separated along the other axis
      intersects(segbbox, boxes[i]) || continue

      edge = edges[i]
      intersection(seg, edge) do I
        itype = type(I)
        # no intersection
        if itype == NotIntersecting
          return
        elseif itype == Overlapping
          # Overlapping segments
          overlap = get(I)
          c, d = vertices(overlap)
          λc = ustrip((c - a) ⋅ v / v²)
          λd = ustrip((d - a) ⋅ v / v²)
          λ₁, λ₂ = minmax(λc, λd)
          push!(segλs, λ₁, λ₂)
          push!(boundaryoverlaps, (λ₁, λ₂))
        else
          # Crossing, EdgeTouching, CornerTouching, etc.
          p = get(I)
          push!(boundarypoints, p)
          λ = ustrip((p - a) ⋅ v / v²)
          push!(segλs, λ)
        end
      end
    end
  end

  # get unique boundary points
  uniquepoints = typeof(a)[]
  for p in boundarypoints
    any(q -> q ≈ p, uniquepoints) || push!(uniquepoints, p)
  end
  boundarypoints = uniquepoints

  # remove duplicate λ values and sort them
  sort!(segλs)
  for i in 2:length(segλs)
    segλs[i] = mayberound(segλs[i], segλs[i - 1])
  end
  unique!(segλs)

  # check if a λ value lies in a boundary overlap
  function inboundaryoverlap(λ)
    any(boundaryoverlaps) do (λ₁, λ₂)
      λ₁ < λ < λ₂ || λ ≈ λ₁ || λ ≈ λ₂
    end
  end

  # create segments from the λ values
  pieces = typeof(seg)[]
  interiorpieces = typeof(seg)[]
  boundarypieces = typeof(seg)[]
  for (λ₁, λ₂) in zip(segλs[1:(end - 1)], segλs[2:end])
    λ₁ ≈ λ₂ && continue
    λmid = (λ₁ + λ₂) / 2
    piece = Segment(seg(λ₁), seg(λ₂))

    if inboundaryoverlap(λmid)
      push!(pieces, piece)
      push!(boundarypieces, piece)
    elseif seg(λmid) ∈ poly
      push!(pieces, piece)
      push!(interiorpieces, piece)
    end
  end

  # merge continuous segments
  function mergepieces(pieces)
    length(pieces) ≤ 1 && return pieces

    merged = eltype(pieces)[]
    for piece in pieces
      if isempty(merged)
        push!(merged, piece)
      else
        previous = last(merged)
        p₁, p₂ = vertices(previous)
        q₁, q₂ = vertices(piece)
        if p₂ ≈ q₁
          merged[end] = Segment(p₁, q₂)
        else
          push!(merged, piece)
        end
      end
    end
    merged
  end

  pieces = mergepieces(pieces)
  interiorpieces = mergepieces(interiorpieces)
  boundarypieces = mergepieces(boundarypieces)

  # create geometry from pieces
  piecegeometry(pieces) = length(pieces) == 1 ? only(pieces) : Multi(pieces)

  # classifying intersections
  if !isempty(interiorpieces)
    geometry = piecegeometry(pieces)
    return @IT Intersecting geometry f # Case 1
  elseif !isempty(boundarypieces)
    geometry = piecegeometry(boundarypieces)
    return @IT EdgeTouching geometry f # Case 2
  elseif length(boundarypoints) == 1
    p = only(boundarypoints)
    if p ≈ a || p ≈ b
      return @IT Touching p f # Case 4
    else
      return @IT CornerTouching p f # Case 3
    end
  elseif !isempty(boundarypoints)
    return @IT Intersecting Multi(boundarypoints) f # Case 1
  else
    return @IT NotIntersecting nothing f # Case 6
  end
end

# Algorithm 4 of Jiménez, J., Segura, R. and Feito, F. 2009.
# (https://www.sciencedirect.com/science/article/pii/S0925772109001448?via%3Dihub)
function intersection(f, seg::Segment{𝔼{3}}, tri::Triangle{𝔼{3}})
  Q1, Q2 = vertices(seg)
  V1, V2, V3 = vertices(tri)

  # according to theorem 1, the algorithm only works
  # when Q1 is not coplanar with the triangle, we need
  # to swap Q1 with Q2 in that case
  if iscoplanar(Q1, V1, V2, V3)
    (Q1, Q2) = (Q2, Q1)
  end

  A = Q1 - V3
  B = V1 - V3
  C = V2 - V3

  W₁ = B × C
  w = A ⋅ W₁

  D = Q2 - V3
  s = D ⋅ W₁

  if w > atol(w)
    # rejection 2
    if s > atol(s)
      return @IT NotIntersecting nothing f
    end

    W₂ = A × D
    t = W₂ ⋅ C

    # rejection 3
    if t < -atol(t)
      return @IT NotIntersecting nothing f
    end

    u = -(W₂ ⋅ B)

    # rejection 4
    if u < -atol(u)
      return @IT NotIntersecting nothing f
    end

    # rejection 5
    if w < (s + t + u)
      return @IT NotIntersecting nothing f
    end
  elseif w < -atol(w)
    # rejection 2
    if s < -atol(s)
      return @IT NotIntersecting nothing f
    end

    W₂ = A × D
    t = W₂ ⋅ C

    # rejection 3
    if t > atol(t)
      return @IT NotIntersecting nothing f
    end

    u = -(W₂ ⋅ B)

    # rejection 4
    if u > atol(u)
      return @IT NotIntersecting nothing f
    end

    # rejection 5
    if w > (s + t + u)
      return @IT NotIntersecting nothing f
    end
  else # w ≈ 0
    if s > atol(s)
      W₂ = D × A
      t = W₂ ⋅ C

      # rejection 3
      if t < -atol(t)
        return @IT NotIntersecting nothing f
      end

      u = -(W₂ ⋅ B)

      # rejection 4
      if u < -atol(u)
        return @IT NotIntersecting nothing f
      end

      # rejection 5
      if -s < (t + u)
        return @IT NotIntersecting nothing f
      end
    elseif s < -atol(s)
      W₂ = D × A
      t = W₂ ⋅ C

      # rejection 3
      if t > atol(t)
        return @IT NotIntersecting nothing f
      end

      u = -(W₂ ⋅ B)

      # rejection 4
      if u > atol(u)
        return @IT NotIntersecting nothing f
      end

      # rejection 5
      if -s > (t + u)
        return @IT NotIntersecting nothing f
      end
    else # s ≈ 0
      # rejection 1, coplanar segment
      return @IT NotIntersecting nothing f
    end
  end

  λ = w / (w - s)
  λ = clamp(λ, zero(λ), one(λ))

  p = Segment(Q1, Q2)(λ)

  return @IT Intersecting p f
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
