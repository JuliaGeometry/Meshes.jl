# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    geodesicfwd(p, ϕ, l)

Solve the direct geodesic problem: return the point reached from the
point `p` after walking the length `l` along the geodesic that leaves
`p` with the azimuth `ϕ`, measured clockwise from the north.

The geodesic is the one of the ellipsoid attached to the datum of `p`,
and is computed with the series of Karney (2013).

See also [`geodesicbwd`](@ref).

## Examples

```julia
p = Point(LatLon(0, 0))

geodesicfwd(p, 90u"°", 1000u"km")

geodesicfwd(p, 90, 1000000)
```

## References

* Karney, C. F. F. 2013. [Algorithms for geodesics](https://doi.org/10.1007/s00190-012-0578-z)
"""
function geodesicfwd(p::Point{🌐}, ϕ, l)
  c = convert(manifoldcrs(p), coords(p))

  # the ellipsoid comes from the datum of the coordinates
  🌎 = ellipsoid(datum(c))

  # the series of Karney need double precision to reach round-off
  T = numtype(lentype(c))
  S = promote_type(T, Float64)
  lat, lon = S(ustrip(c.lat)), S(ustrip(c.lon))
  azi = S(ustrip(u"°", asdeg(ϕ)))
  len = S(ustrip(unit(majoraxis(🌎)), aslen(l)))

  lat′, lon′, _ = _geodesicdirect(🌎, lat, lon, azi, len)

  withcrs(p, (T(lat′), T(lon′)))
end

"""
    geodesicbwd(p₁, p₂)

Solve the inverse geodesic problem: return the azimuth at the point `p₁`
of the shortest geodesic connecting `p₁` to the point `p₂`, measured
clockwise from the north.

The length of that geodesic is [`GeodesicDistance`](@ref), and the
azimuth at the other end is `geodesicbwd(p₂, p₁)`.

See also [`geodesicfwd`](@ref).

## Examples

```julia
p₁ = Point(LatLon(0, 0))
p₂ = Point(LatLon(0, 1))

geodesicbwd(p₁, p₂)
```

## References

* Karney, C. F. F. 2013. [Algorithms for geodesics](https://doi.org/10.1007/s00190-012-0578-z)
"""
function geodesicbwd(p₁::Point{🌐}, p₂::Point{🌐})
  # convert coordinates to same LatLon CRS
  q₁, q₂ = promote(p₁, p₂)
  c₁ = convert(manifoldcrs(q₁), coords(q₁))
  c₂ = convert(manifoldcrs(q₂), coords(q₂))

  🌎 = ellipsoid(datum(c₁))

  # the series of Karney need double precision to reach round-off
  T = numtype(lentype(c₁))
  S = promote_type(T, Float64)
  lat₁, lon₁ = S(ustrip(c₁.lat)), S(ustrip(c₁.lon))
  lat₂, lon₂ = S(ustrip(c₂.lat)), S(ustrip(c₂.lon))

  _, ϕ₁, _ = _geodesicinverse(🌎, lat₁, lon₁, lat₂, lon₂)

  T(ϕ₁) * u"°"
end

"""
    geodesictangent(p, ϕ)

Unit vector tangent to the ellipsoid at the point `p`, pointing along
the azimuth `ϕ`, measured clockwise from the north.

The vector is expressed in the geocentric Cartesian coordinates of the
datum of `p`, and is the direction in which [`geodesicfwd`](@ref) walks.

See also [`geodesicazimuth`](@ref).

## Examples

```julia
p = Point(LatLon(0, 0))

geodesictangent(p, 90u"°")

geodesictangent(p, 90)
```
"""
function geodesictangent(p::Point{🌐}, ϕ)
  T = numtype(lentype(p))
  ê, n̂ = _eastnorth(p)
  s, c = sincosd(T(ustrip(u"°", asdeg(ϕ))))
  unormalize(c * n̂ + s * ê)
end

"""
    geodesicazimuth(p, v)

Azimuth of the vector `v` at the point `p` of the ellipsoid, measured
clockwise from the north. Any component of `v` along the normal of the
ellipsoid is ignored, and the azimuth is undefined if nothing is left.

See also [`geodesictangent`](@ref).

## Examples

```julia
p = Point(LatLon(0, 0))

geodesicazimuth(p, Vec(0, 1, 0))
```
"""
function geodesicazimuth(p::Point{🌐}, v::Vec{3})
  T = numtype(lentype(p))
  ê, n̂ = _eastnorth(p)
  T(atand(v ⋅ ê, v ⋅ n̂)) * u"°"
end

# Solution of the direct geodesic problem: the point reached from (lat₁, lon₁)
# after walking s₁₂ along the geodesic that leaves with azimuth azi₁, following
# Karney (2013) section 3. Angles are in degrees, the length in the unit of the
# semi-axes. Only the position and the azimuth at the second point are returned.
function _geodesicdirect(🌎::Type{<:RevolutionEllipsoid}, lat₁::T, lon₁::T, azi₁::T, s₁₂::T) where {T}
  𝐠 = _geodesicparams(🌎, T)
  tiny = sqrt(floatmin(T))

  # the geodesic is fixed by the azimuth at the first point
  sα₁, cα₁ = sincosd(_anground(_angnormalize(azi₁)))

  # reduced latitude, forced away from the pole so that the quadrants below
  # are unambiguous
  sβ₁, cβ₁ = sincosd(_anground(lat₁))
  sβ₁ *= 𝐠.f₁
  sβ₁, cβ₁ = _hnorm(sβ₁, cβ₁)
  cβ₁ = max(tiny, cβ₁)

  # Clairaut's relation at the northward crossing of the equator, where σ = ω = 0
  sα₀ = sα₁ * cβ₁
  cα₀ = hypot(cα₁, sα₁ * sβ₁)

  sσ₁ = sβ₁
  sω₁ = sα₀ * sβ₁
  cσ₁ = cω₁ = (iszero(sβ₁) && iszero(cα₁)) ? one(T) : cβ₁ * cα₁
  sσ₁, cσ₁ = _hnorm(sσ₁, cσ₁)

  k² = cα₀^2 * 𝐠.e′²
  ε = k² / (2 * (1 + sqrt(1 + k²)) + k²)

  A₁ = 1 + _A1m1(ε)
  C₁ = _C1(ε)
  B₁₁ = _sinseries(sσ₁, cσ₁, C₁)
  sB₁₁, cB₁₁ = sincos(B₁₁)
  sτ₁ = sσ₁ * cB₁₁ + cσ₁ * sB₁₁
  cτ₁ = cσ₁ * cB₁₁ - sσ₁ * sB₁₁

  # the reverted series gives the arc length directly from the distance,
  # with no iteration (Karney, eq. 20 and 21)
  τ₁₂ = s₁₂ / (𝐠.b * A₁)
  sτ₁₂, cτ₁₂ = sincos(τ₁₂)
  B₁₂ = -_sinseries(sτ₁ * cτ₁₂ + cτ₁ * sτ₁₂, cτ₁ * cτ₁₂ - sτ₁ * sτ₁₂, _C1p(ε))
  σ₁₂ = τ₁₂ - (B₁₂ - B₁₁)
  sσ₁₂, cσ₁₂ = sincos(σ₁₂)

  if abs(𝐠.f) > T(0.01)
    # the reverted series is inaccurate for strongly flattened ellipsoids,
    # one Newton step on σ₁₂ recovers the accuracy of the direct series
    sσ₂ = sσ₁ * cσ₁₂ + cσ₁ * sσ₁₂
    cσ₂ = cσ₁ * cσ₁₂ - sσ₁ * sσ₁₂
    δ = A₁ * (σ₁₂ + (_sinseries(sσ₂, cσ₂, C₁) - B₁₁)) - s₁₂ / 𝐠.b
    σ₁₂ -= δ / sqrt(1 + k² * sσ₂^2)
    sσ₁₂, cσ₁₂ = sincos(σ₁₂)
  end

  sσ₂ = sσ₁ * cσ₁₂ + cσ₁ * sσ₁₂
  cσ₂ = cσ₁ * cσ₁₂ - sσ₁ * sσ₁₂

  sβ₂ = cα₀ * sσ₂
  cβ₂ = hypot(sα₀, cα₀ * cσ₂)
  # break the degeneracy of a meridian passing over the pole
  iszero(cβ₂) && ((cβ₂, cσ₂) = (tiny, tiny))
  sα₂, cα₂ = sα₀, cα₀ * cσ₂

  # longitude on the auxiliary sphere, mapped back to the ellipsoid (Karney, eq. 8)
  sω₂, cω₂ = sα₀ * sσ₂, cσ₂
  ω₁₂ = atan(sω₂ * cω₁ - cω₂ * sω₁, cω₂ * cω₁ + sω₂ * sω₁)
  C₃ = _C3(𝐠.c₃, ε)
  B₃₁₂ = _sinseries(sσ₂, cσ₂, C₃) - _sinseries(sσ₁, cσ₁, C₃)
  λ₁₂ = ω₁₂ - 𝐠.f * sα₀ * evalpoly(ε, 𝐠.a₃) * (σ₁₂ + B₃₁₂)

  lat₂ = atand(sβ₂, 𝐠.f₁ * cβ₂)
  lon₂ = _angnormalize(_angnormalize(lon₁) + _angnormalize(rad2deg(λ₁₂)))
  azi₂ = atand(sα₂, cα₂)

  (lat₂, lon₂, azi₂)
end

# Solution of the inverse geodesic problem on an ellipsoid of revolution
# following Karney, C. F. F. 2013. Algorithms for geodesics. Journal of
# Geodesy 87(1), 43-55. https://doi.org/10.1007/s00190-012-0578-z
#
# The series are truncated at order 6, which reaches round-off in double
# precision for ellipsoids with terrestrial flattening. Newton's method is
# used on the azimuth at the first point, with the bracketing interval
# bisected whenever a Newton step leaves the legal range. This keeps the
# solution accurate for nearly antipodal points, where the simpler series
# of Vincenty fails to converge.
#
# Only oblate ellipsoids are handled, which is what CoordRefSystems.jl can
# represent, so the prolate branches of the reference algorithm are omitted.
#
# Return length of the shortest geodesic between two points given in degrees,
# together with the forward azimuths at both points in degrees.
function _geodesicinverse(🌎::Type{<:RevolutionEllipsoid}, lat₁::T, lon₁::T, lat₂::T, lon₂::T) where {T}
  𝐠 = _geodesicparams(🌎, T)
  tiny = sqrt(floatmin(T))
  tol = eps(T)
  maxit₁ = 20
  maxit₂ = maxit₁ + precision(T) + 10

  # longitude difference, made positive
  lon₁₂, δlon₁₂ = _angdiff(lon₁, lon₂)
  lonsign = signbit(lon₁₂) ? -1 : 1
  lon₁₂ *= lonsign
  δlon₁₂ *= lonsign
  λ₁₂ = deg2rad(lon₁₂)
  sλ₁₂, cλ₁₂ = _sincosde(lon₁₂, δlon₁₂)
  lon₁₂ₛ = (T(180) - lon₁₂) - δlon₁₂

  # canonical form 0 ≤ lon₁₂ ≤ 180, -90 ≤ lat₁ ≤ -0 and lat₁ ≤ lat₂ ≤ -lat₁,
  # which the branches below rely on. The three transforms lonsign, swapp and
  # latsign record the reflections and are undone on the azimuths at the end.
  lat₁ = _anground(lat₁)
  lat₂ = _anground(lat₂)
  swapp = abs(lat₁) < abs(lat₂) ? -1 : 1
  if swapp < 0
    lonsign = -lonsign
    lat₁, lat₂ = lat₂, lat₁
  end
  latsign = signbit(lat₁) ? 1 : -1
  lat₁ *= latsign
  lat₂ *= latsign

  # reduced latitudes
  sβ₁, cβ₁ = sincosd(lat₁)
  sβ₁ *= 𝐠.f₁
  sβ₁, cβ₁ = _hnorm(sβ₁, cβ₁)
  cβ₁ = max(tiny, cβ₁)
  sβ₂, cβ₂ = sincosd(lat₂)
  sβ₂ *= 𝐠.f₁
  sβ₂, cβ₂ = _hnorm(sβ₂, cβ₂)
  cβ₂ = max(tiny, cβ₂)

  # force β₂ = ±β₁ exactly when the difference vanishes
  if cβ₁ < -sβ₁
    cβ₂ == cβ₁ && (sβ₂ = copysign(sβ₁, sβ₂))
  else
    abs(sβ₂) == -sβ₁ && (cβ₂ = cβ₁)
  end

  dn₁ = sqrt(1 + 𝐠.e′² * sβ₁^2)
  dn₂ = sqrt(1 + 𝐠.e′² * sβ₂^2)

  s₁₂ = zero(T)
  sα₁ = cα₁ = sα₂ = cα₂ = zero(T)

  # the geodesic may lie on a meridian
  meridian = lat₁ == -90 || iszero(sλ₁₂)
  if meridian
    # head to the target longitude, and arrive heading north
    sα₁, cα₁ = sλ₁₂, cλ₁₂
    sα₂, cα₂ = zero(T), one(T)
    sσ₁, cσ₁ = sβ₁, cα₁ * cβ₁
    sσ₂, cσ₂ = sβ₂, cα₂ * cβ₂
    σ₁₂ = atan(max(zero(T), cσ₁ * sσ₂ - sσ₁ * cσ₂), cσ₁ * cσ₂ + sσ₁ * sσ₂)
    s₁₂, m₁₂, _ = _lengths(𝐠, 𝐠.n, σ₁₂, sσ₁, cσ₁, dn₁, sσ₂, cσ₂, dn₂)
    # zero length geodesics might yield a negative reduced length
    (σ₁₂ < 3tiny || (σ₁₂ < tol && (s₁₂ < 0 || m₁₂ < 0))) && (s₁₂ = zero(T))
    s₁₂ *= 𝐠.b
  end

  if !meridian && iszero(sβ₁) && (iszero(𝐠.f) || lon₁₂ₛ ≥ 𝐠.f * 180)
    # the geodesic runs along the equator
    sα₁ = sα₂ = one(T)
    cα₁ = cα₂ = zero(T)
    s₁₂ = 𝐠.a * λ₁₂
  elseif !meridian
    σ₁₂, sα₁, cα₁, sα₂, cα₂, dnₘ = _inversestart(𝐠, sβ₁, cβ₁, dn₁, sβ₂, cβ₂, dn₂, λ₁₂, sλ₁₂, cλ₁₂)
    if σ₁₂ ≥ 0
      # short line, the auxiliary sphere solution is accurate enough
      s₁₂ = σ₁₂ * 𝐠.b * dnₘ
    else
      # Newton's method on the azimuth α₁, with a bracketing interval that is
      # bisected whenever a Newton step would leave the legal range
      sσ₁ = cσ₁ = sσ₂ = cσ₂ = ε = zero(T)
      sα₁ₐ, cα₁ₐ = tiny, one(T)
      sα₁ᵦ, cα₁ᵦ = tiny, -one(T)
      numit = 0
      tripn = false
      tripb = false
      while true
        v, dv, sα₂, cα₂, σ₁₂, sσ₁, cσ₁, sσ₂, cσ₂, ε =
          _lambda12(𝐠, sβ₁, cβ₁, dn₁, sβ₂, cβ₂, dn₂, sα₁, cα₁, sλ₁₂, cλ₁₂, numit < maxit₁)
        (tripb || !(abs(v) ≥ (tripn ? 8 : 1) * tol) || numit == maxit₂) && break

        # shrink the bracketing interval
        if v > 0 && (numit > maxit₁ || cα₁ / sα₁ > cα₁ᵦ / sα₁ᵦ)
          sα₁ᵦ, cα₁ᵦ = sα₁, cα₁
        elseif v < 0 && (numit > maxit₁ || cα₁ / sα₁ < cα₁ₐ / sα₁ₐ)
          sα₁ₐ, cα₁ₐ = sα₁, cα₁
        end

        newton = false
        if numit < maxit₁ && dv > 0
          dα₁ = -v / dv
          if abs(dα₁) < π
            sdα₁, cdα₁ = sincos(dα₁)
            sα₁′ = sα₁ * cdα₁ + cα₁ * sdα₁
            if sα₁′ > 0
              cα₁ = cα₁ * cdα₁ - sα₁ * sdα₁
              sα₁ = sα₁′
              sα₁, cα₁ = _hnorm(sα₁, cα₁)
              # the slope vanishes in some regimes, so the convergence
              # criterion is based on eps and not on its square root
              tripn = abs(v) ≤ 16tol
              newton = true
            end
          end
        end

        if !newton
          sα₁ = (sα₁ₐ + sα₁ᵦ) / 2
          cα₁ = (cα₁ₐ + cα₁ᵦ) / 2
          sα₁, cα₁ = _hnorm(sα₁, cα₁)
          tripn = false
          tripb = abs(sα₁ₐ - sα₁) + (cα₁ₐ - cα₁) < tol || abs(sα₁ - sα₁ᵦ) + (cα₁ - cα₁ᵦ) < tol
        end

        numit += 1
      end
      s₁₂, _, _ = _lengths(𝐠, ε, σ₁₂, sσ₁, cσ₁, dn₁, sσ₂, cσ₂, dn₂)
      s₁₂ *= 𝐠.b
    end
  end

  # undo the canonicalising transforms on the azimuths
  if swapp < 0
    sα₁, sα₂ = sα₂, sα₁
    cα₁, cα₂ = cα₂, cα₁
  end
  sα₁ *= swapp * lonsign
  cα₁ *= swapp * latsign
  sα₂ *= swapp * lonsign
  cα₂ *= swapp * latsign

  (s₁₂, atand(sα₁, cα₁), atand(sα₂, cα₂))
end

# ------------------------
# MISCELLANEOUS UTILITIES
# ------------------------

# East and north unit vectors of the local frame at the point p, expressed in
# geocentric Cartesian coordinates. Together with the normal of the ellipsoid
# they form the frame in which azimuths are measured.
function _eastnorth(p::Point{🌐})
  c = convert(manifoldcrs(p), coords(p))
  sφ, cφ = sincosd(ustrip(c.lat))
  sλ, cλ = sincosd(ustrip(c.lon))
  u = unit(lentype(p))
  ê = Vec(-sλ * u, cλ * u, zero(sλ) * u)
  n̂ = Vec(-sφ * cλ * u, -sφ * sλ * u, cφ * u)
  (ê, n̂)
end

# quantities that Karney's series need on top of the parameters that
# CoordRefSystems.jl already provides for the ellipsoid of revolution
function _geodesicparams(🌎::Type{<:RevolutionEllipsoid}, ::Type{T}) where {T}
  a = T(ustrip(majoraxis(🌎)))
  b = T(ustrip(minoraxis(🌎)))
  f = T(flattening(🌎))
  e² = T(eccentricity²(🌎))
  f₁ = 1 - f
  e′² = e² / (1 - e²)  # second eccentricity squared
  n = f / (2 - f)      # third flattening
  # threshold on σ₁₂ below which the auxiliary sphere solution is accurate enough
  τ = T(0.1) * sqrt(eps(T)) / sqrt(max(T(0.001), abs(f)) * min(one(T), 1 - f / 2) / 2)
  (; a, b, f, f₁, e′², n, a₃=_A3coeffs(n), c₃=_C3coeffs(n), τ)
end

# distance and reduced length divided by the polar semi-axis (Karney, eq. 15 and 38)
function _lengths(𝐠, ε, σ₁₂, sσ₁, cσ₁, dn₁, sσ₂, cσ₂, dn₂)
  A₁ = _A1m1(ε)
  A₂ = _A2m1(ε)
  C₁ = _C1(ε)
  C₂ = _C2(ε)
  m₀ = A₁ - A₂
  A₁ += 1
  A₂ += 1
  B₁ = _sinseries(sσ₂, cσ₂, C₁) - _sinseries(sσ₁, cσ₁, C₁)
  B₂ = _sinseries(sσ₂, cσ₂, C₂) - _sinseries(sσ₁, cσ₁, C₂)
  J₁₂ = m₀ * σ₁₂ + (A₁ * B₁ - A₂ * B₂)
  s₁₂ = A₁ * (σ₁₂ + B₁)
  m₁₂ = dn₂ * (cσ₁ * sσ₂) - dn₁ * (sσ₁ * cσ₂) - cσ₁ * cσ₂ * J₁₂
  (s₁₂, m₁₂, m₀)
end

# positive root of k⁴ + 2k³ - (x² + y² - 1)k² - 2y²k - y² = 0 (Karney, eq. 65)
function _astroid(x::T, y::T) where {T}
  p = x * x
  q = y * y
  r = (p + q - 1) / 6
  (iszero(q) && r ≤ 0) && return zero(T)
  S = p * q / 4
  r² = r * r
  r³ = r * r²
  disc = S * (S + 2r³)
  u = r
  if disc ≥ 0
    t³ = S + r³
    t³ += t³ < 0 ? -sqrt(disc) : sqrt(disc)
    t = cbrt(t³)
    u += t + (iszero(t) ? zero(T) : r² / t)
  else
    # the cube root is complex but u is real, this form avoids cancellation
    ang = atan(sqrt(-disc), -(S + r³))
    u += 2r * cos(ang / 3)
  end
  v = sqrt(u * u + q)
  uv = u < 0 ? q / (v - u) : u + v
  w = (uv - q) / (2v)
  uv / (sqrt(uv + w * w) + w)
end

# starting azimuth for Newton's method (Karney, section 5)
function _inversestart(𝐠, sβ₁::T, cβ₁, dn₁, sβ₂, cβ₂, dn₂, λ₁₂, sλ₁₂, cλ₁₂) where {T}
  σ₁₂ = -one(T)
  dnₘ = one(T)
  sβ₁₂ = sβ₂ * cβ₁ - cβ₂ * sβ₁
  cβ₁₂ = cβ₂ * cβ₁ + sβ₂ * sβ₁
  sβ₁₂ₐ = sβ₂ * cβ₁ + cβ₂ * sβ₁

  short = cβ₁₂ ≥ 0 && sβ₁₂ < T(0.5) && cβ₂ * λ₁₂ < T(0.5)
  if short
    sβₘ² = (sβ₁ + sβ₂)^2
    sβₘ² /= sβₘ² + (cβ₁ + cβ₂)^2
    dnₘ = sqrt(1 + 𝐠.e′² * sβₘ²)
    sω₁₂, cω₁₂ = sincos(λ₁₂ / (𝐠.f₁ * dnₘ))
  else
    sω₁₂, cω₁₂ = sλ₁₂, cλ₁₂
  end

  sα₁ = cβ₂ * sω₁₂
  cα₁ = cω₁₂ ≥ 0 ? sβ₁₂ + cβ₂ * sβ₁ * sω₁₂^2 / (1 + cω₁₂) : sβ₁₂ₐ - cβ₂ * sβ₁ * sω₁₂^2 / (1 - cω₁₂)

  sσ₁₂ = hypot(sα₁, cα₁)
  cσ₁₂ = sβ₁ * sβ₂ + cβ₁ * cβ₂ * cω₁₂

  sα₂ = zero(T)
  cα₂ = zero(T)
  if short && sσ₁₂ < 𝐠.τ
    # really short line, no iteration needed
    sα₂ = cβ₁ * sω₁₂
    cα₂ = sβ₁₂ - cβ₁ * sβ₂ * (cω₁₂ ≥ 0 ? sω₁₂^2 / (1 + cω₁₂) : 1 - cω₁₂)
    sα₂, cα₂ = _hnorm(sα₂, cα₂)
    σ₁₂ = atan(sσ₁₂, cσ₁₂)
  elseif !(abs(𝐠.n) > T(0.1) || cσ₁₂ ≥ 0 || sσ₁₂ ≥ 6 * abs(𝐠.n) * T(π) * cβ₁^2)
    # nearly antipodal, the spherical guess is useless and the astroid is solved
    # in a frame where the antipodal point is at the origin (Karney, section 5)
    λ₁₂ₓ = atan(-sλ₁₂, -cλ₁₂)
    k² = sβ₁^2 * 𝐠.e′²
    ε = k² / (2 * (1 + sqrt(1 + k²)) + k²)
    λscale = 𝐠.f * cβ₁ * evalpoly(ε, 𝐠.a₃) * T(π)
    βscale = λscale * cβ₁
    x = λ₁₂ₓ / λscale
    y = sβ₁₂ₐ / βscale
    if y > -200eps(T) && x > -1 - 1000sqrt(eps(T))
      # strip near the cut
      sα₁ = min(one(T), -x)
      cα₁ = -sqrt(1 - sα₁^2)
    else
      k = _astroid(x, y)
      ω₁₂ₐ = -λscale * x * k / (1 + k)
      sω₁₂ = sin(ω₁₂ₐ)
      cω₁₂ = -cos(ω₁₂ₐ)
      sα₁ = cβ₂ * sω₁₂
      cα₁ = sβ₁₂ₐ - cβ₂ * sβ₁ * sω₁₂^2 / (1 - cω₁₂)
    end
  end

  if sα₁ > 0
    sα₁, cα₁ = _hnorm(sα₁, cα₁)
  else
    sα₁, cα₁ = one(T), zero(T)
  end

  (σ₁₂, sα₁, cα₁, sα₂, cα₂, dnₘ)
end

# residual λ₁₂(α₁) - λ₁₂ and its derivative (Karney, section 4)
function _lambda12(𝐠, sβ₁::T, cβ₁, dn₁, sβ₂, cβ₂, dn₂, sα₁, cα₁, sλ₁₂, cλ₁₂, diff) where {T}
  # break the degeneracy of the equatorial line
  (iszero(sβ₁) && iszero(cα₁)) && (cα₁ = -sqrt(floatmin(T)))

  # Clairaut's relation sinα₀ = sinα cosβ
  sα₀ = sα₁ * cβ₁
  cα₀ = hypot(cα₁, sα₁ * sβ₁)

  sσ₁ = sβ₁
  sω₁ = sα₀ * sβ₁
  cσ₁ = cω₁ = cα₁ * cβ₁
  sσ₁, cσ₁ = _hnorm(sσ₁, cσ₁)

  sα₂ = cβ₂ ≠ cβ₁ ? sα₀ / cβ₂ : sα₁
  cα₂ = if cβ₂ ≠ cβ₁ || abs(sβ₂) ≠ -sβ₁
    sqrt((cα₁ * cβ₁)^2 + (cβ₁ < -sβ₁ ? (cβ₂ - cβ₁) * (cβ₁ + cβ₂) : (sβ₁ - sβ₂) * (sβ₁ + sβ₂))) / cβ₂
  else
    abs(cα₁)
  end

  sσ₂ = sβ₂
  sω₂ = sα₀ * sβ₂
  cσ₂ = cω₂ = cα₂ * cβ₂
  sσ₂, cσ₂ = _hnorm(sσ₂, cσ₂)

  σ₁₂ = atan(max(zero(T), cσ₁ * sσ₂ - sσ₁ * cσ₂), cσ₁ * cσ₂ + sσ₁ * sσ₂)
  sω₁₂ = max(zero(T), cω₁ * sω₂ - sω₁ * cω₂)
  cω₁₂ = cω₁ * cω₂ + sω₁ * sω₂
  η = atan(sω₁₂ * cλ₁₂ - cω₁₂ * sλ₁₂, cω₁₂ * cλ₁₂ + sω₁₂ * sλ₁₂)

  k² = cα₀^2 * 𝐠.e′²
  ε = k² / (2 * (1 + sqrt(1 + k²)) + k²)
  C₃ = _C3(𝐠.c₃, ε)
  B₃ = _sinseries(sσ₂, cσ₂, C₃) - _sinseries(sσ₁, cσ₁, C₃)
  v = η - 𝐠.f * evalpoly(ε, 𝐠.a₃) * sα₀ * (σ₁₂ + B₃)

  # the derivative follows from the reduced length (Karney, eq. 38)
  dv = if !diff
    zero(T)
  elseif iszero(cα₂)
    -2 * 𝐠.f₁ * dn₁ / sβ₁
  else
    _, m₁₂, _ = _lengths(𝐠, ε, σ₁₂, sσ₁, cσ₁, dn₁, sσ₂, cσ₂, dn₂)
    m₁₂ * 𝐠.f₁ / (cα₂ * cβ₂)
  end

  (v, dv, sα₂, cα₂, σ₁₂, sσ₁, cσ₁, sσ₂, cσ₂, ε)
end

# ----------------
# ANGLE UTILITIES
# ----------------

_hnorm(x, y) = (h = hypot(x, y); (x / h, y / h))

# coarsen x to the resolution available near 1/16 degree, so that angles
# that are tiny but not zero do not turn into near singular cases. The
# expression is not a no-op, the rounding happens in the two subtractions.
function _anground(x::T) where {T}
  z = T(1) / 16
  y = abs(x)
  w = z - y
  y = w > 0 ? z - w : y
  copysign(y, x)
end

# reduce x to (-180, 180]
function _angnormalize(x::T) where {T}
  y = rem(x, T(360), RoundNearest)
  abs(y) == 180 ? copysign(T(180), x) : y
end

# sum with error term
function _twosum(u::T, v::T) where {T}
  s = u + v
  uᵣ = s - v
  vᵣ = s - uᵣ
  uᵣ -= u
  vᵣ -= v
  (s, iszero(s) ? s : -(uᵣ + vᵣ))
end

# difference y - x in (-180, 180] with error term
function _angdiff(x::T, y::T) where {T}
  d, e = _twosum(rem(-x, T(360), RoundNearest), rem(y, T(360), RoundNearest))
  d, e = _twosum(rem(d, T(360), RoundNearest), e)
  if iszero(d) || abs(d) == 180
    d = copysign(d, iszero(e) ? y - x : -e)
  end
  (d, e)
end

# sine and cosine of (x + t) degrees with exact reduction of x
function _sincosde(x::T, t::T) where {T}
  d = rem(x, T(90), RoundNearest)
  q = round(Int, (x - d) / 90)
  d = _anground(d + t)
  r = deg2rad(d)
  s, c = sincos(r)
  if 2abs(d) == 90
    c = sqrt(T(1) / 2)
    s = copysign(c, r)
  elseif 3abs(d) == 90
    c = sqrt(T(3)) / 2
    s = copysign(T(1) / 2, r)
  end
  m = q & 3
  m == 0 ? (s, c) : m == 1 ? (c, -s) : m == 2 ? (-s, -c) : (-c, s)
end

# Clenshaw summation of ∑ cₗ sin(2lx)
@inline function _sinseries(sinx::T, cosx::T, c::NTuple{N,T}) where {N,T}
  α = 2 * (cosx - sinx) * (cosx + sinx)
  b₁ = zero(T)
  b₂ = zero(T)
  for cₗ in reverse(c)
    b₁, b₂ = α * b₁ - b₂ + cₗ, b₁
  end
  2 * sinx * cosx * b₁
end

# --------------
# KARNEY SERIES
# --------------

# A₁ - 1 (Karney, eq. 17)
_A1m1(ε::T) where {T} = (evalpoly(ε * ε, (T(0), T(64), T(4), T(1))) / 256 + ε) / (1 - ε)

# A₂ - 1 (Karney, eq. 42)
_A2m1(ε::T) where {T} = (evalpoly(ε * ε, (T(0), T(-192), T(-28), T(-11))) / 256 - ε) / (1 + ε)

# C₁ₗ (Karney, eq. 18)
function _C1(ε::T) where {T}
  ε² = ε * ε
  (
    ε * evalpoly(ε², (T(-16), T(6), T(-1))) / 32,
    ε^2 * evalpoly(ε², (T(-128), T(64), T(-9))) / 2048,
    ε^3 * evalpoly(ε², (T(-16), T(9))) / 768,
    ε^4 * evalpoly(ε², (T(-5), T(3))) / 512,
    ε^5 * T(-7) / 1280,
    ε^6 * T(-7) / 2048
  )
end

# C₁ₗ′, the reverted series taking τ back to σ (Karney, eq. 21)
function _C1p(ε::T) where {T}
  ε² = ε * ε
  (
    ε * evalpoly(ε², (T(768), T(-432), T(205))) / 1536,
    ε^2 * evalpoly(ε², (T(3840), T(-4736), T(4005))) / 12288,
    ε^3 * evalpoly(ε², (T(116), T(-225))) / 384,
    ε^4 * evalpoly(ε², (T(2695), T(-7173))) / 7680,
    ε^5 * T(3467) / 7680,
    ε^6 * T(38081) / 61440
  )
end

# C₂ₗ (Karney, eq. 43)
function _C2(ε::T) where {T}
  ε² = ε * ε
  (
    ε * evalpoly(ε², (T(16), T(2), T(1))) / 32,
    ε^2 * evalpoly(ε², (T(384), T(64), T(35))) / 2048,
    ε^3 * evalpoly(ε², (T(80), T(15))) / 768,
    ε^4 * evalpoly(ε², (T(35), T(7))) / 512,
    ε^5 * T(63) / 1280,
    ε^6 * T(77) / 2048
  )
end

# coefficients of A₃ as a polynomial in ε (Karney, eq. 24)
function _A3coeffs(n::T) where {T}
  (
    one(T),
    (n - 1) / 2,
    evalpoly(n, (T(-2), T(-1), T(3))) / 8,
    -evalpoly(n, (T(1), T(3), T(1))) / 16,
    -(2n + 3) / 64,
    T(-3) / 128
  )
end

# coefficients of C₃ₗ as polynomials in ε (Karney, eq. 25)
function _C3coeffs(n::T) where {T}
  (
    (
      evalpoly(n, (T(1), T(-1))) / 4,
      evalpoly(n, (T(1), T(0), T(-1))) / 8,
      evalpoly(n, (T(3), T(3), T(-1))) / 64,
      evalpoly(n, (T(5), T(2))) / 128,
      T(3) / 128
    ),
    (
      evalpoly(n, (T(2), T(-3), T(1))) / 32,
      evalpoly(n, (T(3), T(-2), T(-3))) / 64,
      evalpoly(n, (T(3), T(1))) / 128,
      T(5) / 256
    ),
    (evalpoly(n, (T(5), T(-9), T(5))) / 192, evalpoly(n, (T(9), T(-10))) / 384, T(7) / 512),
    (evalpoly(n, (T(7), T(-14))) / 512, T(7) / 512),
    (T(21) / 2560,)
  )
end

function _C3(c₃, ε::T) where {T}
  (
    ε * evalpoly(ε, c₃[1]),
    ε^2 * evalpoly(ε, c₃[2]),
    ε^3 * evalpoly(ε, c₃[3]),
    ε^4 * evalpoly(ε, c₃[4]),
    ε^5 * evalpoly(ε, c₃[5])
  )
end
