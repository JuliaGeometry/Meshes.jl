# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ∠(A, B, C)

Angle ∠ABC between rays BA and BC.
See https://en.wikipedia.org/wiki/Angle.

Vec2 angles return a value in range [-π, π].
Vec3 angles return a value in range [0, π].

## Examples

```julia
∠(Point(1,0), Point(0,0), Point(0,1)) == π/2
```
"""
∠(A::P, B::P, C::P) where {P<:Point{2}} = ∠(A-B, C-B)
∠(A::P, B::P, C::P) where {P<:Point{3}} = ∠(A-B, C-B)

"""
    ∠(u, v)

Angle between vectors u and v.
See https://en.wikipedia.org/wiki/Angle.

Vec2 angles return a value in range [-π, π].
Vec3 angles return a value in range [0, π].

## Examples

```julia
∠(Vec(1,0), Vec(0,1)) == π/2
```

Thank you Jeffrey Sarnoff for contributing this work.
"""
function ∠(u::V, v::V) where {V<:Vec}
    T = eltype(u.coords)

    u_unit = unitize(u)
    v_unit = unitize(v)

    y = norm2(u_unit .- v_unit)
    x = norm2(u_unit .+ v_unit)

    a = 2 * atan(y, x)

    a = !(signbit(a) || signbit(pi - a)) ? a : (signbit(a) ? zero(T) : (T)(pi))

    ifelse(isnegangle(u_unit, v_unit), -a, a)
end

@inline norm2(u::V) where {V<:Vec} = sqrt(foldl(+, abs2.(u.coords)))
@inline unitize(u::V) where {V<:Vec} = u.coords ./ norm2(u)

@inline isnegangle(u::V, v::V) where {V<:Vec3} = false
@inline isnegangle(u::V, v::V) where {V<:Vec2} =
    u.coords[1] * v.coords[2] < u.coords[2] * v.coords[1]

