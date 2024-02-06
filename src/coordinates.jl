# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Coordinates{N}

Parent type of all coordinate types.
"""
abstract type Coordinates{N} end

Base.isapprox(c₁::C, c₂::C; kwargs...) where {C<:Coordinates} =
  all(isapprox(getfield(c₁, n), getfield(c₂, n); kwargs...) for n in fieldnames(C))

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, coords::Coordinates)
  name = nameof(typeof(coords))
  print(io, "$name(")
  printfields(io, coords, compact=true)
  print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", coords::Coordinates)
  name = nameof(typeof(coords))
  print(io, "$name coordinates")
  printfields(io, coords)
end

# ----------
# UTILITIES
# ----------

const Len{T} = Quantity{T,u"𝐋"}
const Rad{T} = Quantity{T,NoDims,typeof(u"rad")}
const Deg{T} = Quantity{T,NoDims,typeof(u"°")}

# only add the unit if the argument is not a quantity
addunit(x::Number, u) = x * u
addunit(x::Quantity, u) = throw(ArgumentError("invalid units for coordinates, please check the documentation"))

# ----------------
# IMPLEMENTATIONS
# ----------------

include("coordinates/basic.jl")
include("coordinates/gis.jl")

# ------------
# CONVERSIONS
# ------------

# Cartesian <-> Polar
Base.convert(::Type{<:Cartesian}, (; ρ, ϕ)::Polar) = Cartesian(ρ * cos(ϕ), ρ * sin(ϕ))
function Base.convert(::Type{<:Polar}, (; coords)::Cartesian{2})
  x, y = coords
  Polar(sqrt(x^2 + y^2), atanpos(y, x) * u"rad")
end

# Cartesian <-> Cylindrical
Base.convert(::Type{<:Cartesian}, (; ρ, ϕ, z)::Cylindrical) = Cartesian(ρ * cos(ϕ), ρ * sin(ϕ), z)
function Base.convert(::Type{<:Cylindrical}, (; coords)::Cartesian{3})
  x, y, z = coords
  Cylindrical(sqrt(x^2 + y^2), atanpos(y, x) * u"rad", z)
end

# Cartesian <-> Spherical
Base.convert(::Type{<:Cartesian}, (; r, θ, ϕ)::Spherical) =
  Cartesian(r * sin(θ) * cos(ϕ), r * sin(θ) * sin(ϕ), r * cos(θ))
function Base.convert(::Type{<:Spherical}, (; coords)::Cartesian{3})
  x, y, z = coords
  Spherical(sqrt(x^2 + y^2 + z^2), atan(sqrt(x^2 + y^2), z) * u"rad", atanpos(y, x) * u"rad")
end

# WGS84 ellipsoid
const a = 6378137.0
const f⁻¹ = 298.257223563
const f = inv(f⁻¹)
const b = a * (1 - f)
const e² = (2 - f) / f⁻¹
const e = √e²

# LatLon <-> Mercator
function fwdformula(::Type{Mercator})
  x = (ϕ, λ) -> oftype(λ, a) * λ
  y = (ϕ, λ) -> begin
    ey = oftype(ϕ, e)
    oftype(ϕ, a) * (asinh(tan(ϕ)) - ey * atanh(ey * sin(ϕ)))
  end
  x, y
end

function Base.convert(::Type{Mercator}, (; coords)::LatLon)
  lat, lon = coords
  ϕ = ustrip(deg2rad(lat))
  λ = ustrip(deg2rad(lon))
  x, y = fwdformula(Mercator)
  Mercator(x(ϕ, λ) * u"m", y(ϕ, λ) * u"m")
end

# LatLon <-> WebMercator
fwdformula(::Type{WebMercator}) = ((ϕ, λ) -> oftype(λ, a) * λ, (ϕ, λ) -> oftype(ϕ, a) * asinh(tan(ϕ)))
invformula(::Type{WebMercator}) = ((x, y) -> atan(sinh(y / oftype(y, a))), (x, y) -> x / oftype(x, a))

function Base.convert(::Type{WebMercator}, (; coords)::LatLon)
  lat, lon = coords
  ϕ = ustrip(deg2rad(lat))
  λ = ustrip(deg2rad(lon))
  x, y = fwdformula(WebMercator)
  WebMercator(x(ϕ, λ) * u"m", y(ϕ, λ) * u"m")
end

function Base.convert(::Type{LatLon}, (; coords)::WebMercator)
  x, y = coords
  nx = ustrip(x)
  ny = ustrip(y)
  ϕ, λ = invformula(WebMercator)
  LatLon(rad2deg(ϕ(nx, ny)) * u"°", rad2deg(λ(nx, ny)) * u"°")
end

# adjust negative angles
function atanpos(y, x)
  α = atan(y, x)
  ifelse(α ≥ zero(α), α, α + oftype(α, 2π))
end
