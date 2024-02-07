# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# ----------
# UTILITIES
# ----------

const Len{T} = Quantity{T,u"𝐋"}
const Met{T} = Quantity{T,u"𝐋",typeof(u"m")}
const Rad{T} = Quantity{T,NoDims,typeof(u"rad")}
const Deg{T} = Quantity{T,NoDims,typeof(u"°")}

# only add the unit if the argument is not a quantity
addunit(x::Number, u) = x * u
addunit(x::Quantity, u) = throw(ArgumentError("invalid units for coordinates, please check the documentation"))

# adjust negative angles
function atanpos(y, x)
  α = atan(y, x)
  ifelse(α ≥ zero(α), α, α + oftype(α, 2π))
end

# -----------
# ELLIPSOIDS
# -----------

const WGS84 = let
  a = 6378137.0 * u"m"
  f⁻¹ = 298.257223563
  f = inv(f⁻¹)
  b = a * (1 - f)
  e² = (2 - f) / f⁻¹
  e = √e²
  (; a, b, e, e², f, f⁻¹)
end

# ------------
# COORDINATES
# ------------

"""
    Coordinates{N}

Parent type of all coordinate types.
"""
abstract type Coordinates{N} end

_fields(coords::Coordinates) = coords
_fnames(coords::Coordinates) = fieldnames(typeof(coords))

function Base.isapprox(coords₁::C, coords₂::C; kwargs...) where {N,C<:Coordinates{N}}
  f₁ = _fields(coords₁)
  f₂ = _fields(coords₂)
  all(ntuple(i -> isapprox(getfield(f₁, i), getfield(f₂, i); kwargs...), N))
end

# -----------
# IO METHODS
# -----------

function Base.show(io::IO, coords::Coordinates)
  name = prettyname(coords)
  print(io, "$name(")
  fields = _fields(coords)
  fnames = _fnames(coords)
  printfields(io, fields, fnames, compact=true)
  print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", coords::Coordinates)
  name = prettyname(coords)
  print(io, "$name coordinates")
  fields = _fields(coords)
  fnames = _fnames(coords)
  printfields(io, fields, fnames)
end

#------
# EPSG
#------

"""
    EPSG{code}

EPSG dataset `code` between 1024 and 32767.
Codes can be searched at [epsg.io](https://epsg.io/).

See [EPSG Geodetic Parameter Dataset](https://en.wikipedia.org/wiki/EPSG_Geodetic_Parameter_Dataset)
"""
struct EPSG{Code,N,Coords} <: Coordinates{N}
  coords::Coords
end

EPSG{Code,N,Coords}(args...) where {Code,N,Coords} = EPSG{Code,N,Coords}(Coords(args))

_fields(coords::EPSG) = coords.coords
_fnames(coords::EPSG) = keys(coords.coords)

"""
    typealias(::Type{EPSG{code}})

Returns a coordinate type that has the EPSG `code`.
"""
function typealias end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("coordinates/basic.jl")
include("coordinates/latlon.jl")
include("coordinates/mercator.jl")
include("coordinates/webmercator.jl")
include("coordinates/platecarree.jl")

# ----------
# FALLBACKS
# ----------

Base.convert(T::Type{EPSG{Code}}, coords::Coordinates) where {Code} = convert(typealias(T), coords)
