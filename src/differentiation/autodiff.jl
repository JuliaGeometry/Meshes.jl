# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function deriv(geom::Geometry, uvw, j, method::DI.AbstractADType)
  dim = paramdim(geom)
  function f(t)
    args = ntuple(i -> i == j ? t : uvw[i], dim)
    p = geom(args...)
    CoordRefSystems.raw(coords(p))
  end
  dvals = DI.derivative(f, method.backend, uvw[j])
  ℒ = lentype(geom)
  u = unit(ℒ)
  Vec(Tuple(dvals .* u))
end
