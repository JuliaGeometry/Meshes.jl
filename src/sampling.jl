# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Sampler

A method for sampling from geometries.
"""
abstract type Sampler end

"""
    sample(geometry, sampler)

Sample elements from `geometry` with `sampler`.
"""
sample(geometry::Geometry, sizes) =
  sample(geometry, RegularSampler(sizes))

include("sampling/regular.jl")
