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
  m₁, M₁ = to.(extrema(box₁))
  m₂, M₂ = to.(extrema(box₂))

  # relevant vertices
  u = withcrs(box₁, max.(promote(m₁, m₂)...))
  v = withcrs(box₁, min.(promote(M₁, M₂)...))

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
  m₁, M₁ = coords.(extrema(box₁ |> Proj(crs)))
  m₂, M₂ = coords.(extrema(box₂ |> Proj(crs)))

  # check latitude coordinate
  latₛ, latₑ = max(m₁.lat, m₂.lat), min(M₁.lat, M₂.lat)
  latₛ > latₑ && return @IT NotIntersecting nothing f

  # check longitude coordinate
  lons = _lonintersects(m₁.lon, M₁.lon, m₂.lon, M₂.lon)
  isempty(lons) && return @IT NotIntersecting nothing f

  # classify
  latpoint = isapproxzero(latₑ - latₛ)
  lonpoints = all(r -> isapproxzero(last(r) - first(r)), lons)
  corner = latpoint && lonpoints
  overlap = !latpoint && !lonpoints

  # construct geometry
  geoms = map(lons) do (lonₛ, lonₑ)
    u = withcrs(box₁, (latₛ, lonₛ))
    if corner
      u
    else
      v = withcrs(box₁, (latₑ, lonₑ))
      Box(u, v)
    end
  end
  geom = length(geoms) == 1 ? only(geoms) : Multi(geoms)

  if overlap
    return @IT Overlapping geom f
  elseif corner
    return @IT CornerTouching geom f
  else
    return @IT Touching geom f
  end
end

function _lonintersects(lo₁, hi₁, lo₂, hi₂)
  ranges₁ = _lonranges(lo₁, hi₁)
  ranges₂ = _lonranges(lo₂, hi₂)
  ranges = Tuple{typeof(lo₁),typeof(lo₁)}[]
  for (a, b) in ranges₁, (c, d) in ranges₂
    lo, hi = max(a, c), min(b, d)
    lo ≤ hi && push!(ranges, (lo, hi))
  end

  # -180° and 180° represent the same meridian
  Δ = oftype(lo₁, 180u"°")
  if _containsantimeridian(ranges₁, Δ) && _containsantimeridian(ranges₂, Δ) && !_containsantimeridian(ranges, Δ)
    push!(ranges, (Δ, Δ))
  end

  # pieces meeting at the antimeridian form a single wrapped interval
  if length(ranges) ≥ 2 && first(ranges)[1] == -Δ && last(ranges)[2] == Δ
    lo, hi = last(ranges)[1], first(ranges)[2]
    return lo == Δ && hi == -Δ ? [(Δ, Δ)] : [(lo, hi)]
  end
  ranges
end

# longitude intervals in the canonical [-180°, 180°] range
function _lonranges(lo, hi)
  Δ = oftype(lo, 180u"°")
  lo ≤ hi ? ((lo, hi),) : ((-Δ, hi), (lo, Δ))
end

_containsantimeridian(ranges, Δ) = any(r -> first(r) == -Δ || last(r) == Δ, ranges)
