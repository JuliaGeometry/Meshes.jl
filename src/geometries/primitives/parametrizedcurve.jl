# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ParametrizedCurve(func, range = (0.0, 1.0))

A parametrized curve is a curve defined by a function `func` that maps a parameter `t` to a `Point` in space.
The parameter `t` is defined in the interval `[a, b]`. The curve can only be evaluated for `t` in the
interval `[0, 1]`.

## Examples

```julia
ParametrizedCurve(t -> Point(cospi(2t), sinpi(2t)))
```
"""
struct ParametrizedCurve{M<:Meshes.Manifold,C<:Meshes.CRS,T<:Real,F<:Function} <: Primitive{M,C}
  a::T
  b::T
  func::F

  function ParametrizedCurve(func, ab=(0.0, 1.0))
    a, b = promote(ab...)
    p = func(a)
    M = manifold(p)
    C = crs(p)
    T = typeof(a)
    new{M,C,T,typeof(func)}(a, b, func)
  end
end

paramdim(::ParametrizedCurve) = 1
startparameter(curve::ParametrizedCurve) = curve.a
Base.minimum(curve::ParametrizedCurve) = curve(0.0)
endparameter(curve::ParametrizedCurve) = curve.b
Base.maximum(curve::ParametrizedCurve) = curve(1.0)
Base.extrema(curve::ParametrizedCurve) = curve(0.0), curve(1.0)

function (curve::ParametrizedCurve)(t)
  if t < 0.0 || t > 1.0
    throw(DomainError(t, "c(t) is not defined for t outside [0, 1]."))
  end
  a, b = curve.range
  curve.func(a + t * (b - a))
end
