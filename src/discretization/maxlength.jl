# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    MaxLengthDiscretization(length)

TODO
"""
struct MaxLengthDiscretization{ℒ<:Len} <: DiscretizationMethod
  length::ℒ
  MaxLengthDiscretization(length::ℒ) where {ℒ<:Len} = new{float(ℒ)}(length)
end

MaxLengthDiscretization(length) = MaxLengthDiscretization(addunit(length, u"m"))

function discretize(box::Box, method::MaxLengthDiscretization)
  sizes = ceil.(Int, sides(box) ./ method.length)
  discretize(box, RegularDiscretization(sizes))
end

function discretize(segment::Segment, method::MaxLengthDiscretization)
  size = ceil(Int, _measure(segment) / method.length)
  discretize(segment, RegularDiscretization(size))
end

discretize(chain::Chain, method::MaxLengthDiscretization) =
  mapreduce(s -> discretize(s, method), merge, segments(chain))

# -----------------
# HELPER FUNCTIONS
# -----------------

_measure(s::Segment{<:𝔼}) = measure(s)

# TODO: Haversine returns the shortest distance between two points
# this is not always equal to the distance between two directed points
function _measure(s::Segment{<:🌐})
  T = numtype(lentype(s))
  🌎 = ellipsoid(datum(crs(s)))
  r = numconvert(T, majoraxis(🌎))

  a, b = extrema(s)
  evaluate(Haversine(r), a, b)
end
