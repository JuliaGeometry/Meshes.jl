# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

abstract type Datum end

abstract type NoDatum <: Datum end

ellipsoid(::Type{NoDatum}) = nothing
latitudeₒ(::Type{NoDatum}) = nothing
longitudeₒ(::Type{NoDatum}) = nothing
altitudeₒ(::Type{NoDatum}) = nothing

abstract type WGS84 <: Datum end

ellipsoid(::Type{WGS84}) = WGS84🌎
latitudeₒ(::Type{WGS84}) = 0.0u"°"
longitudeₒ(::Type{WGS84}) = 0.0u"°"
altitudeₒ(::Type{WGS84}) = 0.0u"m"
