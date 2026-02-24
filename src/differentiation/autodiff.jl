# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    AutoDiff(backend=DI.AutoMooncakeForward())

Automatic differentiation method using DifferentiationInterface.jl.

The default backend is MooncakeForward. Any DI-compatible backend can be used:

    AutoDiff()                              # MooncakeForward (default)
    AutoDiff(DI.AutoMooncakeForward())      # MooncakeForward (explicit)
"""
struct AutoDiff{B} <: DifferentiationMethod
  backend::B
end

AutoDiff() = AutoDiff(DI.AutoMooncakeForward())

function deriv(geom::Geometry, uvw, j, method::AutoDiff)
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
