# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MaxLengthDiscretization(length)

Discretize geometries into parts with sides of maximum `length` in length units (default to meters).
"""
struct MaxLengthDiscretization{ℒ<:Len} <: DiscretizationMethod
  length::ℒ
  MaxLengthDiscretization(length::ℒ) where {ℒ<:Len} = new{float(ℒ)}(length)
end

MaxLengthDiscretization(length) = MaxLengthDiscretization(aslen(length))

_discretize(geometry::Geometry, method::MaxLengthDiscretization) =
  refine(discretize(geometry), MaxLengthRefinement(method.length))

function _discretize(box::Box, method::MaxLengthDiscretization)
  sizes = ceil.(Int, approxsides(box) ./ method.length)
  discretize(box, RegularDiscretization(sizes))
end

function _discretize(segment::Segment, method::MaxLengthDiscretization)
  size = ceil(Int, measure(segment) / method.length)
  discretize(segment, RegularDiscretization(size))
end

_discretize(chain::Chain, method::MaxLengthDiscretization) =
  mapreduce(s -> _discretize(s, method), merge, segments(chain))

function _discretize(quad::Quadrangle, method::MaxLengthDiscretization)
  sizes = ceil.(Int, approxsides(quad) ./ method.length)
  discretize(quad, RegularDiscretization(sizes))
end

function _discretize(hexa::Hexahedron, method::MaxLengthDiscretization)
  sizes = ceil.(Int, approxsides(hexa) ./ method.length)
  discretize(hexa, RegularDiscretization(sizes))
end
