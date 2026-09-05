# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# The intersection type can be one of four types:
# 1. overlap with non-zero measure (Overlapping -> Box)
# 2. intersect at corner point (CornerTouching -> Point)
# 3. intersect at one of the facets (Touching -> Box)
# 4. do not overlap nor intersect (NotIntersecting -> Nothing)
function intersection(f, box₁::Box, box₂::Box)
  # retrieve corner points
  min₁, max₁ = to.(extrema(box₁))
  min₂, max₂ = to.(extrema(box₂))

  # relevant vertices
  u = withcrs(box₁, max.(promote(min₁, min₂)...))
  v = withcrs(box₁, min.(promote(max₁, max₂)...))

  # auxiliary variables
  δ = v - u
  δ̄ = abs.(δ)
  τ = atol(eltype(δ))

  # branch on possible configurations
  if all(>(τ), δ)
    return @IT Overlapping Box(u, v) f
  elseif all(<(τ), δ̄)
    return @IT CornerTouching u f
  elseif any(<(τ), δ̄) && (δ == δ̄ || δ == -δ̄)
    return @IT Touching (u ⪯ v ? Box(u, v) : Box(v, u)) f
  else
    return @IT NotIntersecting nothing f
  end
end

# The intersection type can be one of four types:
# 1. overlap with non-zero measure (Overlapping -> Box)
# 2. intersect at corner point (CornerTouching -> Point)
# 3. intersect at one of the facets (Touching -> Box)
# 4. do not overlap nor intersect (NotIntersecting -> Nothing)
function intersection(f, box₁::Box{🌐}, box₂::Box{🌐})
  # corners in a common CRS
  crs = manifoldcrs(box₁)
  min₁, max₁ = coords.(extrema(box₁ |> Proj(crs)))
  min₂, max₂ = coords.(extrema(box₂ |> Proj(crs)))

  # check latitude coordinate
  latₛ, latₑ = max(min₁.lat, min₂.lat), min(max₁.lat, max₂.lat)
  latₛ > latₑ && return @IT NotIntersecting nothing f

  # check longitude coordinate
  lonranges = _lonintersects(min₁.lon, max₁.lon, min₂.lon, max₂.lon)
  isempty(lonranges) && return @IT NotIntersecting nothing f

  # classify intersection
  singlelat = isapproxequal(latₛ, latₑ)
  singlelon = all(lon -> isapproxequal(lon[1], lon[2]), lonranges)
  iscorner = singlelat && singlelon
  overlaps = !singlelat && !singlelon

  # construct geometry
  geoms = map(lonranges) do (lonₛ, lonₑ)
    u = withcrs(box₁, (latₛ, lonₛ))
    v = withcrs(box₁, (latₑ, lonₑ))
    iscorner ? u : Box(u, v)
  end
  geom = length(geoms) == 1 ? only(geoms) : Multi(geoms)

  if overlaps
    return @IT Overlapping geom f
  elseif iscorner
    return @IT CornerTouching geom f
  else
    return @IT Touching geom f
  end
end

function _lonintersects(minlon₁::T, maxlon₁::T, minlon₂::T, maxlon₂::T) where {T}
  # -180° and 180° represent the same meridian
  πdeg = T(180u"°")

  # convert range to one or two non-decreasing ranges,
  # depending on whether it crosses the antimeridian
  ranges₁ = minlon₁ ≤ maxlon₁ ? ((minlon₁, maxlon₁),) : ((-πdeg, maxlon₁), (minlon₁, πdeg))
  ranges₂ = minlon₂ ≤ maxlon₂ ? ((minlon₂, maxlon₂),) : ((-πdeg, maxlon₂), (minlon₂, πdeg))

  # compute the intersection of the two sets of ranges
  ranges = Tuple{T,T}[]
  for (a₁, b₁) in ranges₁, (a₂, b₂) in ranges₂
    a, b = max(a₁, a₂), min(b₁, b₂)
    a ≤ b && push!(ranges, (a, b))
  end

  # check if the intersection contains the antimeridian
  _hasπdeg(ranges) = any(lon -> lon[1] == -πdeg || lon[2] == πdeg, ranges)
  if _hasπdeg(ranges₁) && _hasπdeg(ranges₂) && !_hasπdeg(ranges)
    push!(ranges, (πdeg, πdeg))
  end

  # ranges meeting at the antimeridian form a single wrapped interval
  if length(ranges) ≥ 2 && first(ranges)[1] == -πdeg && last(ranges)[2] == πdeg
    a, b = last(ranges)[1], first(ranges)[2]
    return a == πdeg && b == -πdeg ? [(πdeg, πdeg)] : [(a, b)]
  end

  ranges
end
