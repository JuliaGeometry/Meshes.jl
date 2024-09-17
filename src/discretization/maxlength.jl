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
  sizes = ceil.(Int, _sides(box) ./ method.length)
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

_sides(box::Box{<:𝔼}) = sides(box)

function _sides(box::Box{<:🌐})
  A, B = extrema(box)
  a = convert(LatLon, coords(A))
  b = convert(LatLon, coords(B))
  P = withcrs(box, (a.lat, b.lon), LatLon)

  AP = Segment(A, P)
  PB = Segment(P, B)
  (_measure(AP), _measure(PB))
end

_measure(segment::Segment{<:𝔼}) = measure(segment)

# TODO: Haversine returns the shortest distance between two points
# this is not always equal to the distance between two directed points
function _measure(segment::Segment{<:🌐})
  T = numtype(lentype(segment))
  🌎 = ellipsoid(datum(crs(segment)))
  r = numconvert(T, majoraxis(🌎))

  A, B = extrema(segment)
  evaluate(Haversine(r), A, B)
end
