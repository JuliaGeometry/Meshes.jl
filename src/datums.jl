# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

abstract type Datum end

abstract type NoDatum <: Datum end

ellipsoid(::Type{NoDatum}) = nothing
latitude₀(::Type{NoDatum}) = nothing
longitude₀(::Type{NoDatum}) = nothing
altitude₀(::Type{NoDatum}) = nothing

abstract type WGS84 <: Datum end

ellipsoid(::Type{WGS84}) = WGS84🌎
latitude₀(::Type{WGS84}) = 0.0u"°"
longitude₀(::Type{WGS84}) = 0.0u"°"
altitude₀(::Type{WGS84}) = 0.0u"m"
