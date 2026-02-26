# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ParametrizedCurve(fun, range = (0.0, 1.0))

A parametrized curve is a curve defined by a function `fun` that maps a
(unitless) parameter `t` in the given `range` to a `Point` in space.

## Examples

```julia
ParametrizedCurve(t -> Point(cos(t), sin(t)), (0, 2π))
```
"""
struct ParametrizedCurve{M<:Manifold,C<:CRS,F<:Function,R<:Tuple} <: Primitive{M,C}
  fun::F
  range::R
  ParametrizedCurve{M,C}(fun::F, range::R) where {M<:Manifold,C<:CRS,F<:Function,R<:Tuple} = new{M,C,F,R}(fun, range)
end

function ParametrizedCurve(fun, range=(0.0, 1.0))
  a, b = promote(range...)
  r = (a, b)
  p = fun(a)
  ParametrizedCurve{manifold(p),crs(p)}(fun, r)
end

paramdim(::Type{<:ParametrizedCurve}) = 1

Base.minimum(curve::ParametrizedCurve) = curve.fun(first(curve.range))

Base.maximum(curve::ParametrizedCurve) = curve.fun(last(curve.range))

Base.extrema(curve::ParametrizedCurve) = minimum(curve), maximum(curve)

function (curve::ParametrizedCurve)(t)
  a, b = curve.range
  curve.fun(a + t * (b - a))
end
