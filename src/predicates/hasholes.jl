# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    hasholes(geometry)

Tells whether or not the `geometry` contains holes.
"""
hasholes(g::Geometry) = hasholes(typeof(g))

hasholes(::Type{<:Geometry}) = false

hasholes(p::Pgon) = length(rings(p)) > 1
