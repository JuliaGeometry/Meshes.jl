# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GeometricDistance

A distance between points in geometries.
"""
abstract type GeometricDistance end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("distances/euclidean.jl")
include("distances/mahalanobis.jl")
include("distances/geodesic.jl")
