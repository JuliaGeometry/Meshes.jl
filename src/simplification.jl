# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SimplificationMethod

A method for simplifying geometries.
"""
abstract type SimplificationMethod end

"""
    simplify(geometry, method)

Simplify `geometry` with given `method`.
"""
function simplify end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("simplification/douglaspeucker.jl")
