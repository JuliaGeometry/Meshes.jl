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
struct ParametrizedCurve{M<:Manifold,C<:CRS,F<:Function,R<:Tuple} <: Primitive{M,C}
  func::F
  range::R

  function ParametrizedCurve(func, range=(0.0, 1.0))
    a = first(range)
    p = func(a)
    M = manifold(p)
    C = crs(p)
    new{M,C,typeof(func),typeof(range)}(func, range)
  end
end

paramdim(::ParametrizedCurve) = 1
Base.minimum(curve::ParametrizedCurve) = curve(0.0)
Base.maximum(curve::ParametrizedCurve) = curve(1.0)
Base.extrema(curve::ParametrizedCurve) = curve(0.0), curve(1.0)

function (curve::ParametrizedCurve)(t)
  if t < 0.0 || t > 1.0
    throw(DomainError(t, "c(t) is not defined for t outside [0, 1]."))
  end
  a, b = curve.range
  curve.func(a + t * (b - a))
end
