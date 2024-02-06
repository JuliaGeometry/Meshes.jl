# Cartesian <-> Polar
Base.convert(::Type{Cartesian}, (; ρ, ϕ)::Polar) = Cartesian(ρ * cos(ϕ), ρ * sin(ϕ))
function Base.convert(::Type{Polar}, (; coords)::Cartesian{2})
  x, y = coords
  Polar(sqrt(x^2 + y^2), atanpos(y, x) * u"rad")
end

# Cartesian <-> Cylindrical
Base.convert(::Type{Cartesian}, (; ρ, ϕ, z)::Cylindrical) = Cartesian(ρ * cos(ϕ), ρ * sin(ϕ), z)
function Base.convert(::Type{Cylindrical}, (; coords)::Cartesian{3})
  x, y, z = coords
  Cylindrical(sqrt(x^2 + y^2), atanpos(y, x) * u"rad", z)
end

# Cartesian <-> Spherical
Base.convert(::Type{Cartesian}, (; r, θ, ϕ)::Spherical) =
  Cartesian(r * sin(θ) * cos(ϕ), r * sin(θ) * sin(ϕ), r * cos(θ))
function Base.convert(::Type{Spherical}, (; coords)::Cartesian{3})
  x, y, z = coords
  Spherical(sqrt(x^2 + y^2 + z^2), atan(sqrt(x^2 + y^2), z) * u"rad", atanpos(y, x) * u"rad")
end

# LatLon <-> Mercator
function Base.convert(::Type{Mercator}, (; coords)::LatLon)
  λ = ustrip(deg2rad(coords.lon))
  ϕ = ustrip(deg2rad(coords.lat))
  aλ = oftype(λ, wgs84.a)
  aϕ = oftype(ϕ, wgs84.a)
  eϕ = oftype(ϕ, wgs84.e)
  x = aλ * λ
  y = aϕ * (asinh(tan(ϕ)) - eϕ * atanh(eϕ * sin(ϕ)))
  Mercator(x * u"m", y * u"m")
end

# LatLon <-> WebMercator
function Base.convert(::Type{WebMercator}, (; coords)::LatLon)
  λ = ustrip(deg2rad(coords.lon))
  ϕ = ustrip(deg2rad(coords.lat))
  aλ = oftype(λ, wgs84.a)
  aϕ = oftype(ϕ, wgs84.a)
  x = aλ * λ
  y = aϕ * asinh(tan(ϕ))
  WebMercator(x * u"m", y * u"m")
end

function Base.convert(::Type{LatLon}, (; coords)::WebMercator)
  x = ustrip(coords.x)
  y = ustrip(coords.y)
  ax = oftype(x, wgs84.a)
  ay = oftype(y, wgs84.a)
  λ = x / ax
  ϕ = atan(sinh(y / ay))
  LatLon(rad2deg(ϕ) * u"°", rad2deg(λ) * u"°")
end

# adjust negative angles
function atanpos(y, x)
  α = atan(y, x)
  ifelse(α ≥ zero(α), α, α + oftype(α, 2π))
end
