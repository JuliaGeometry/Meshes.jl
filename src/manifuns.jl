# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

check_point(::Type{𝔼{Dim}}, p::Point; kwargs...) where {Dim} = embeddim(p) == Dim

function check_point(::Type{🌐}, p::Point; kwargs...)
  check = iszeroalt(coords(p); kwargs...)
  check ? nothing : DomainError(p, "Point $p is not on the ellispoid 🌐")
end

function iszeroalt(coords::CRS; kwargs...)
  lla = convert(LatLonAlt, coords)
  isapprox(lla.alt, zero(lla.alt); kwargs...)
end
