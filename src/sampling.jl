# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    SamplingMethod

A method for sampling from geometries.
"""
abstract type SamplingMethod end

"""
    sample(geometry, method)

Sample elements from `geometry` with `method`.
"""
sample(geometry::Geometry, sizes) =
  sample(geometry, RegularSampling(sizes))

include("sampling/regular.jl")
