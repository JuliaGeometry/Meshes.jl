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

function discretize(segment::Segment{<:𝔼}, method::MaxLengthDiscretization)
  size = ceil(Int, length(segment) / method.length)
  discretize(segment, RegularDiscretization(size))
end

function discretize(segment::Segment{<:🌐}, method::MaxLengthDiscretization)
  T = numtype(lentype(segment))
  🌎 = ellipsoid(datum(crs(segment)))
  r = numconvert(T, majoraxis(🌎))

  a, b = vertices(segment)
  d = evaluate(Haversine(r), a, b)

  size = ceil(Int, d / method.length)

  discretize(segment, RegularDiscretization(size))
end

discretize(chain::Chain, method::MaxLengthDiscretization) =
  mapreduce(s -> discretize(s, method), merge, segments(chain))
