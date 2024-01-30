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
