# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ClippingMethod

A method for clipping polygons with other polygons.
"""
abstract type ClippingMethod end

"""
    clip(subject, other, method)

Clip the `subject` polygon with `other` polygon using clipping `method`.
"""
function clip end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("clipping/sutherlandhodgman.jl")
include("clipping/greinerhormann.jl")
