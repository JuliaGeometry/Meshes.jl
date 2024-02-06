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

include("coordinates/ellipsoids.jl")
include("coordinates/basic.jl")
include("coordinates/gis.jl")

# ------------
# CONVERSIONS
# ------------

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
  lat, lon = coords
  λ = ustrip(deg2rad(lon))
  ϕ = ustrip(deg2rad(lat))
  x = oftype(λ, a) * λ
  ey = oftype(ϕ, e)
  y = oftype(ϕ, a) * (asinh(tan(ϕ)) - ey * atanh(ey * sin(ϕ)))
  Mercator(x * u"m", y * u"m")
end

# LatLon <-> WebMercator
function Base.convert(::Type{WebMercator}, (; coords)::LatLon)
  lat, lon = coords
  λ = ustrip(deg2rad(lon))
  ϕ = ustrip(deg2rad(lat))
  x = oftype(λ, a) * λ
  y = oftype(ϕ, a) * asinh(tan(ϕ))
  WebMercator(x * u"m", y * u"m")
end

function Base.convert(::Type{LatLon}, (; coords)::WebMercator)
  x, y = coords
  nx = ustrip(x)
  ny = ustrip(y)
  λ = nx / oftype(nx, a)
  ϕ = atan(sinh(ny / oftype(ny, a)))
  LatLon(rad2deg(ϕ) * u"°", rad2deg(λ) * u"°")
end

# adjust negative angles
function atanpos(y, x)
  α = atan(y, x)
  ifelse(α ≥ zero(α), α, α + oftype(α, 2π))
end
