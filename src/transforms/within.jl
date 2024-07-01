# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

struct Within{T} <: GeometricTransform
  bounds::T
end

_addunit(x::Nothing) = x
_addunit(x) = addunit.(x, u"m")

Within(; x=nothing, y=nothing, z=nothing) = Within((_addunit(x), _addunit(y), _addunit(z)))

function _applybound(dim, bound, box)
  Dim = embeddim(box)
  if Dim < dim || isnothing(bound)
    box
  else
    bmin, bmax = bound
    min = to(minimum(box))
    max = to(maximum(box))
    nmin = Vec(ntuple(i -> i == dim ? bmin : min[i], Dim))
    nmax = Vec(ntuple(i -> i == dim ? bmax : max[i], Dim))
    Box(withdatum(box, nmin), withdatum(box, nmax))
  end
end

function preprocess(t::Within, d::Domain)
  bbox = boundingbox(d)
  for (dim, bound) in enumerate(t.bounds)
    bbox = _applybound(dim, bound, bbox)
  end
  indices(d, bbox)
end

function apply(t::Within, d::Domain)
  inds = preprocess(t, d)
  n = view(d, inds)
  n, nothing
end
