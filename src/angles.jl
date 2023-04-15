# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ∠(A, B, C)

Angle ∠ABC between rays BA and BC.
See https://en.wikipedia.org/wiki/Angle.

Return a value in range [-π, π] for 2D
angles and a value in range [0, π] for
3D angles.

## Examples

```julia
∠(Point(1,0), Point(0,0), Point(0,1)) == π/2
```
"""
∠(A::Point, B::Point, C::Point) = ∠(A-B, C-B)

"""
    ∠(u, v)

Angle between vectors u and v.
See https://en.wikipedia.org/wiki/Angle.

Return a value in range [-π, π] for 2D
angles and a value in range [0, π] for
3D angles.

## Examples

```julia
∠(Vec(1,0), Vec(0,1)) == π/2
```
"""
function ∠(u::Vec, v::Vec)
  T = eltype(u)

  uunit = unitize(u)
  vunit = unitize(v)

  y = norm2(uunit .- vunit)
  x = norm2(uunit .+ vunit)

  a = 2 * atan(y, x)

  a = !(signbit(a) || signbit(pi - a)) ? a : (signbit(a) ? zero(T) : (T)(pi))

  ifelse(isnegangle(uunit, vunit), -a, a)
end

@inline norm2(u::Vec) = sqrt(foldl(+, abs2.(u)))

@inline unitize(u::Vec) = u ./ norm2(u)

@inline isnegangle(u::Vec{2}, v::Vec{2}) = u[1] * v[2] < u[2] * v[1]
@inline isnegangle(u::Vec{3}, v::Vec{3}) = false