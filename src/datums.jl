# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Datum

Parent type of all datum types.
"""
abstract type Datum end

"""
    ellipsoid(D)

Returns the ellipsoid of the datum type `D`.
"""
function ellipsoid end

"""
    latitudeₒ(D)

Returns the latitude origin of the datum type `D`.
"""
function latitudeₒ end

"""
    longitudeₒ(D)

Returns the longitude origin of the datum type `D`.
"""
function longitudeₒ end

"""
    altitudeₒ(D)

Returns the altitude origin of the datum type `D`.
"""
function altitudeₒ end

"""
    NoDatum

Represents the absence of datum in a CRS.
"""
abstract type NoDatum <: Datum end

ellipsoid(::Type{NoDatum}) = nothing
latitudeₒ(::Type{NoDatum}) = nothing
longitudeₒ(::Type{NoDatum}) = nothing
altitudeₒ(::Type{NoDatum}) = nothing

"""
    WGS84

WGS84 (World Geodetic System) datum.

See [World Geodetic System](https://en.wikipedia.org/wiki/World_Geodetic_System)
"""
abstract type WGS84 <: Datum end

ellipsoid(::Type{WGS84}) = WGS84🌎
latitudeₒ(::Type{WGS84}) = 0.0u"°"
longitudeₒ(::Type{WGS84}) = 0.0u"°"
altitudeₒ(::Type{WGS84}) = 0.0u"m"

"""
    WGS84

Winkel Tripel datum.

See [ESRI:53042](https://epsg.io/53042)
"""
abstract type WIII <: Datum end

ellipsoid(::Type{WIII}) = WIII🌎
latitudeₒ(::Type{WIII}) = 0.0u"°"
longitudeₒ(::Type{WIII}) = 0.0u"°"
altitudeₒ(::Type{WIII}) = 0.0u"m"
