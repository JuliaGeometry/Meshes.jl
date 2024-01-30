# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Coordinates{N,T}

Parent type of all coordinate types.
"""
abstract type Coordinates{N,T} end

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

# ----------------
# IMPLEMENTATIONS
# ----------------

include("coordinates/basic.jl")
include("coordinates/gis.jl")

# ------------
# CONVERSIONS
# ------------

# Cartesian <--> Polar
Base.convert(::Type{<:Cartesian}, (; ρ, ϕ)::Polar) = Cartesian(ρ * cos(ϕ), ρ * sin(ϕ))
function Base.convert(::Type{<:Polar}, (; coords)::Cartesian{2})
  x, y = coords
  # adjust the atan(y, x) interval: [-π,π] + π = [0,2π]
  Polar(sqrt(x^2 + y^2), atan(y, x) + π)
end

# Cartesian <--> Cylindrical
Base.convert(::Type{<:Cartesian}, (; ρ, ϕ, z)::Cylindrical) = Cartesian(ρ * cos(ϕ), ρ * sin(ϕ), z)
function Base.convert(::Type{<:Cylindrical}, (; coords)::Cartesian{3})
  x, y, z = coords
  # adjust the atan(y, x) interval: [-π,π] + π = [0,2π]
  Cylindrical(sqrt(x^2 + y^2), atan(y, x) + π, z)
end

# Cartesian <--> Spherical
Base.convert(::Type{<:Cartesian}, (; r, θ, ϕ)::Spherical) = Cartesian(r * sin(θ) * cos(ϕ), r * sin(θ) * sin(ϕ), r * cos(θ))
function Base.convert(::Type{<:Spherical}, (; coords)::Cartesian{3})
  x, y, z = coords
  # adjust the atan(y, x) interval: [-π,π] + π = [0,2π]
  Spherical(sqrt(x^2 + y^2 + z^2), atan(sqrt(x^2 + y^2), z), atan(y, x) + π)
end
