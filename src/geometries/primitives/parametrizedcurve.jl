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
ParametrizedCurve(t -> Point(cos(t), sin(t)), (0, 2π))
```
"""
struct ParametrizedCurve{M<:Manifold,C<:CRS,F<:Function,R<:Tuple} <: Primitive{M,C}
  func::F
  range::R

  function ParametrizedCurve(func, range=(0.0, 1.0))
    a, b = promote(range...)
    _range = (a, b)
    p = func(a)
    M = manifold(p)
    C = crs(p)
    new{M,C,typeof(func),typeof(_range)}(func, _range)
  end
end

paramdim(::Type{<:ParametrizedCurve}) = 1
Base.minimum(curve::ParametrizedCurve) = curve.func(first(curve.range))
Base.maximum(curve::ParametrizedCurve) = curve.func(last(curve.range))
Base.extrema(curve::ParametrizedCurve) = minimum(curve), maximum(curve)

function (curve::ParametrizedCurve)(t)
  if t < 0 || t > 1
    throw(DomainError(t, "c(t) is not defined for t outside [0, 1]."))
  end
  a, b = curve.range
  curve.func(a + t * (b - a))
end
