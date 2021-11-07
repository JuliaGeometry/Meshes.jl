# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Neighborhood

A neighborhood is a geometry that is not attached to any specific
point in space, and is free to slide over a domain of interest.
"""
abstract type Neighborhood end

"""
    MetricBall

A metric ball is a neighborhood that can be expressed in terms
of a metric and a set of radii.
"""
abstract type MetricBall <: Neighborhood end

"""
    metric(ball)

Return the metric of the metric `ball`.
"""
function metric end

"""
    radii(ball)

Return the radii of the metric `ball`.
"""
function radii end

# ----------------
# IMPLEMENTATIONS
# ----------------

include("neighborhoods/isotropicball.jl")
include("neighborhoods/anisotropicball.jl")
