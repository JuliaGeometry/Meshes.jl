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

MaxLengthDiscretization(length) = MaxLengthDiscretization(addunit(length, u"m"))

function discretize(box::Box, method::MaxLengthDiscretization)
  sizes = ceil.(Int, _sides(box) ./ method.length)
  discretize(box, RegularDiscretization(sizes))
end

function discretize(segment::Segment, method::MaxLengthDiscretization)
  size = ceil(Int, measure(segment) / method.length)
  discretize(segment, RegularDiscretization(size))
end

discretize(chain::Chain, method::MaxLengthDiscretization) =
  mapreduce(s -> discretize(s, method), merge, segments(chain))

function discretize(geometry::Geometry, method::MaxLengthDiscretization)
  iscoarse(e) = perimeter(e) > method.length * nvertices(e)
  mesh = simplexify(geometry)
  while any(iscoarse, mesh)
    mesh = refine(mesh, TriSubdivision())
  end
  mesh
end

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
