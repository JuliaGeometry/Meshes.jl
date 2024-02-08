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

# ----
# CRS
# ----

struct CRS{ID,Coords,Datum}
  coords::Coords
end

CRS{ID,Coords,Datum}(args...) where {ID,Coords,Datum} = CRS{ID,Coords,Datum}(Coords(args))

_coords(coords::CRS) = getfield(coords, :coords)

Base.propertynames(coords::CRS) = propertynames(_coords(coords))

Base.getproperty(coords::CRS, name::Symbol) = getproperty(_coords(coords), name)

function Base.isapprox(coords₁::C, coords₂::C; kwargs...) where {C<:CRS}
  c₁ = _coords(coords₁)
  c₂ = _coords(coords₂)
  N = length(c₁)
  all(ntuple(i -> isapprox(getfield(c₁, i), getfield(c₂, i); kwargs...), N))
end

# ------
# DATUM
# ------

"""
    datum(coords)

Returns the datum of the coordinates `coords`.
"""
datum(::CRS{ID,Coords,Datum}) where {ID,Coords,Datum} = Datum

"""
    ellipsoid(coords)

Returns the ellipsoid of the coordinates `coords`.
"""
ellipsoid(coords::CRS) = ellipsoid(datum(coords))

"""
    latitudeₒ(coords)

Returns the latitude origin of the coordinates `coords`.
"""
latitudeₒ(coords::CRS) = latitudeₒ(datum(coords))

"""
    longitudeₒ(coords)

Returns the longitude origin of the coordinates `coords`.
"""
longitudeₒ(coords::CRS) = longitudeₒ(datum(coords))

"""
    altitudeₒ(coords)

Returns the altitude origin of the coordinates `coords`.
"""
altitudeₒ(coords::CRS) = altitudeₒ(datum(coords))

# -----------
# IO METHODS
# -----------

_fnames(coords::CRS) = fieldnames(typeof(_coords(coords)))

function Base.show(io::IO, coords::CRS)
  name = prettyname(coords)
  print(io, "$name(")
  printfields(io, _coords(coords), _fnames(coords), compact=true)
  print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", coords::CRS)
  name = prettyname(coords)
  print(io, "$name coordinates")
  printfields(io, _coords(coords), _fnames(coords))
end

#-----------
# EPSG/ESRI
#-----------

"""
    EPSG{code}

EPSG dataset `code` between 1024 and 32767.
Codes can be searched at [epsg.io](https://epsg.io/).

See [EPSG Geodetic Parameter Dataset](https://en.wikipedia.org/wiki/EPSG_Geodetic_Parameter_Dataset)
"""
abstract type EPSG{code} end

"""
    ESRI{code}

ESRI dataset `code`. Codes can be searched at [epsg.io](https://epsg.io/).
"""
abstract type ESRI{code} end

"""
    typealias(::Type{EPSG{code}})
    typealias(::Type{ESRI{code}})

Returns a coordinate type that has the EPSG/ESRI `code`.
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

Base.convert(T::Type{EPSG{Code}}, coords::CRS) where {Code} = convert(typealias(T), coords)
Base.convert(T::Type{ESRI{Code}}, coords::CRS) where {Code} = convert(typealias(T), coords)
