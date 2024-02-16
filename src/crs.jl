# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

#-----------
# EPSG/ESRI
#-----------

"""
    EPSG{code}

EPSG dataset `code` between 1024 and 32767.
Codes can be searched at [epsg.io](https://epsg.io).

See [EPSG Geodetic Parameter Dataset](https://en.wikipedia.org/wiki/EPSG_Geodetic_Parameter_Dataset)
"""
abstract type EPSG{Code} end

"""
    ESRI{code}

ESRI dataset `code`. Codes can be searched at [epsg.io](https://epsg.io).
"""
abstract type ESRI{Code} end

"""
    typealias(::Type{EPSG{code}})
    typealias(::Type{ESRI{code}})

Returns a CRS type that has the EPSG/ESRI `code`.
"""
function typealias end

# ----
# CRS
# ----

"""
    CRS{Datum}

Parent type of all Coordinate Reference Systems (CRS) with a given `Datum`
that can be used to georeference any point in physical space.
"""
abstract type CRS{Datum} end

_coords(coords::CRS) = coords

function Base.isapprox(coords₁::C, coords₂::C; kwargs...) where {C<:CRS}
  c₁ = _coords(coords₁)
  c₂ = _coords(coords₂)
  N = nfields(c₁)
  all(ntuple(i -> isapprox(getfield(c₁, i), getfield(c₂, i); kwargs...), N))
end

# ------
# DATUM
# ------

"""
    datum(coords)

Returns the datum of the coordinates `coords`.
"""
datum(coords::CRS) = datum(typeof(coords))
datum(::Type{<:CRS{Datum}}) where {Datum} = Datum

"""
    ellipsoid(coords)

Returns the ellipsoid of the coordinates `coords`.
"""
ellipsoid(coords::CRS) = ellipsoid(typeof(coords))
ellipsoid(T::Type{<:CRS}) = ellipsoid(datum(T))

"""
    latitudeₒ(coords)

Returns the latitude origin of the coordinates `coords`.
"""
latitudeₒ(coords::CRS) = latitudeₒ(typeof(coords))
latitudeₒ(T::Type{<:CRS}) = latitudeₒ(datum(T))

"""
    longitudeₒ(coords)

Returns the longitude origin of the coordinates `coords`.
"""
longitudeₒ(coords::CRS) = longitudeₒ(typeof(coords))
longitudeₒ(T::Type{<:CRS}) = longitudeₒ(datum(T))

"""
    altitudeₒ(coords)

Returns the altitude origin of the coordinates `coords`.
"""
altitudeₒ(coords::CRS) = altitudeₒ(typeof(coords))
altitudeₒ(T::Type{<:CRS}) = altitudeₒ(datum(T))

# -----------
# IO METHODS
# -----------

_fnames(coords::CRS) = fieldnames(typeof(coords))

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

# ----------------
# IMPLEMENTATIONS
# ----------------

const Len{T} = Quantity{T,u"𝐋"}
const Met{T} = Quantity{T,u"𝐋",typeof(u"m")}
const Rad{T} = Quantity{T,NoDims,typeof(u"rad")}
const Deg{T} = Quantity{T,NoDims,typeof(u"°")}

include("crs/basic.jl")
include("crs/latlon.jl")
include("crs/mercator.jl")
include("crs/webmercator.jl")
include("crs/eqdistcylindrical.jl")
include("crs/eqareacylindrical.jl")
include("crs/winkeltripel.jl")
include("crs/robinson.jl")
include("crs/orthographic.jl")

# ----------
# FALLBACKS
# ----------

Base.convert(T::Type{EPSG{Code}}, coords::CRS) where {Code} = convert(typealias(T), coords)
Base.convert(T::Type{ESRI{Code}}, coords::CRS) where {Code} = convert(typealias(T), coords)
