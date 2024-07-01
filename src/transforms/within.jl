# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

struct Within{X,Y,Z} <: GeometricTransform
  x::X
  y::Y
  z::Z
end

Within(; x=nothing, y=nothing, z=nothing) = Within(
  isnothing(x) ? x : addunit.(x, u"m"),
  isnothing(y) ? y : addunit.(y, u"m"),
  isnothing(z) ? z : addunit.(z, u"m")
)

function preprocess(t::Within, d::Domain)
  bbox = boundingbox(d)
  bbox₁ = _within(1, t.x, bbox)
  bbox₂ = _within(2, t.y, bbox₁)
  bbox₃ = _within(3, t.z, bbox₂)
  indices(d, bbox₃)
end

function apply(t::Within, d::Domain)
  inds = preprocess(t, d)
  n = view(d, inds)
  n, nothing
end

function _within(dim, bound, box)
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
