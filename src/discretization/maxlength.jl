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

function _discretize(box::Box, method::MaxLengthDiscretization)
  sizes = ceil.(Int, _sides(box) ./ method.length)
  discretize(box, RegularDiscretization(sizes))
end

function _discretize(segment::Segment, method::MaxLengthDiscretization)
  size = ceil(Int, measure(segment) / method.length)
  discretize(segment, RegularDiscretization(size))
end

_discretize(chain::Chain, method::MaxLengthDiscretization) =
  mapreduce(s -> _discretize(s, method), merge, segments(chain))

_discretize(geometry::Geometry, method::MaxLengthDiscretization) =
  refine(discretize(geometry), MaxLengthRefinement(method.length))

# -----------------
# HELPER FUNCTIONS
# -----------------

_sides(box::Box{<:𝔼}) = sides(box)

function _sides(box::Box{<:🌐})
  A, B = extrema(box)
  a = convert(LatLon, coords(A))
  b = convert(LatLon, coords(B))
  P = withcrs(box, (a.lat, b.lon), LatLon)
  AP = Segment(A, P)
  PB = Segment(P, B)
  (measure(AP), measure(PB))
end
