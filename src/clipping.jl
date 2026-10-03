# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ClippingMethod

A method for clipping subject geometries with other geometries.
"""
abstract type ClippingMethod end

"""
    clip(subject, other, method)

Clip the `subject` geometry with `other` geometry using clipping `method`.
"""
clip(subject::Geometry, other::Geometry, method::ClippingMethod) = clip(_aspolygon(subject), _aspolygon(other), method)

# ----------------
# IMPLEMENTATIONS
# ----------------

include("clipping/sutherlandhodgman.jl")
include("clipping/greinerhormann.jl")

# -----------------
# HELPER FUNCTIONS
# -----------------

_aspolygon(p::Polygon) = p
_aspolygon(b::Box) = convert(Quadrangle, b)
