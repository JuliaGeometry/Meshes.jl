# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ParametrizedCurve(a, b, func)

A parametrized curve is a curve defined by a function `func` that maps a parameter `t` to a `Point` in space.
The parameter `t` is defined in the interval `[a, b]`.

## Examples

```julia
ParametrizedCurve(0.0, 1.0, t -> Point(cospi(2t), sinpi(2t)))
```
"""
struct ParametrizedCurve{M<:Meshes.Manifold,C<:Meshes.CRS,T<:Real,F<:Function} <: Primitive{M,C}
  a::T
  b::T
  func::F

  function ParametrizedCurve(a::T, b::T, func) where {T<:Real}
    p = func(a)
    M = manifold(p)
    C = crs(p)
    new{M,C,T,typeof(func)}(a, b, func)
  end
end

paramdim(::ParametrizedCurve) = 1
startparameter(curve::ParametrizedCurve) = curve.a
Base.minimum(curve::ParametrizedCurve) = curve(curve.a)
endparameter(curve::ParametrizedCurve) = curve.b
Base.maximum(curve::ParametrizedCurve) = curve(curve.b)
# Base.extrema(curve::ParametrizedCurve) = curve(curve.a), curve(curve.b)

function (curve::ParametrizedCurve)(t)
  if t < startparameter(curve) || t > endparameter(curve)
    throw(DomainError(t, "c(t) is not defined for t outside [$(startparameter(curve)), $(endparameter(curve))]."))
  end
  curve.func(t)
end
