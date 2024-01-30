# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

struct Cartesian{N,T} <: Coordinates{N,T}
  coords::NTuple{N,T}
  Cartesian{N,T}(coords) where {N,T} = new{N,float(T)}(coords)
end

Cartesian(coords::NTuple{N,T}) where {N,T} = Coordinates{N,T}(coords)
Cartesian(coords::Tuple) = Coordinates(promote(coords...))
Cartesian(coords...) = Cartesian(coords)

struct Polar{T} <: Coordinates{2,T}
  radius::T
  angle::T
  Polar{T}(radius, angle) where {T} = new{float(T)}(radius, angle)
end

Polar(radius::T, angle::T) where {T} = Polar{T}(radius, angle)
Polar(radius, angle) = Polar(promote(radius, angle)...)

struct Spherical{T} <: Coordinates{3,T}
  radius::T
  polar::T
  azimuth::T
  Spherical{T}(radius, polar, azimuth) where {T} = new{float(T)}(radius, polar, azimuth)
end

Spherical(radius::T, polar::T, azimuth::T) where {T} = Spherical{T}(radius, polar, azimuth)
Spherical(radius, polar, azimuth) = Spherical(promote(radius, polar, azimuth)...)

struct Cylindrical{T} <: Coordinates{3,T}
  radius::T
  angle::T
  height::T
  Cylindrical{T}(radius, angle, height) where {T} = new{float(T)}(radius, angle, height)
end

Cylindrical(radius::T, angle::T, height::T) where {T} = Cylindrical{T}(radius, angle, height)
Cylindrical(radius, angle, height) = Cylindrical(promote(radius, angle, height)...)
