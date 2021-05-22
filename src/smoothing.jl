# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SmoothingMethod

A method for smoothing meshes.
"""
abstract type SmoothingMethod end

"""
    smooth(mesh, method)

Smooth `mesh` with given method.
"""
function smooth end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("smoothing/taubin.jl")
