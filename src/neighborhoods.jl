# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Neighborhood

A neighborhood is a geometry that is not attached to any specific
point in the space, and is free to slide over a domain of interest.
"""
abstract type Neighborhood end

"""
    MetricBall

A metric ball is a neighborhood that can be expressed in terms of a
metric and a range. They are useful for fast searches with tree
data structures.
"""
abstract type MetricBall <: Neighborhood end

"""
    metric(ball)

Return the metric of the norm `ball`.
"""
function metric(::MetricBall) end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("neighborhoods/normball.jl")
include("neighborhoods/ellipsoid.jl")
