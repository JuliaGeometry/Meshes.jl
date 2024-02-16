# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

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

# ----------------
# IMPLEMENTATIONS
# ----------------

typealias(::Type{EPSG{3395}}) = Mercator{WGS84}
typealias(::Type{EPSG{3857}}) = WebMercator{WGS84}
typealias(::Type{EPSG{4326}}) = LatLon{WGS84}
typealias(::Type{EPSG{32662}}) = PlateCarree{WGS84}
typealias(::Type{ESRI{54017}}) = Behrmann{WGS84}
typealias(::Type{ESRI{54030}}) = Robinson{WGS84}
typealias(::Type{ESRI{54034}}) = Lambert{WGS84}
typealias(::Type{ESRI{54042}}) = WinkelTripel{WGS84}
typealias(::Type{ESRI{102035}}) = Orthographic{90.0u"°",0.0u"°",true,WGS84}
typealias(::Type{ESRI{102037}}) = Orthographic{-90.0u"°",0.0u"°",true,WGS84}

# ----------
# FALLBACKS
# ----------

Base.convert(T::Type{EPSG{Code}}, coords::CRS) where {Code} = convert(typealias(T), coords)
Base.convert(T::Type{ESRI{Code}}, coords::CRS) where {Code} = convert(typealias(T), coords)
