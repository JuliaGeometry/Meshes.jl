# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ∠(A, B, C)

Angle ∠ABC between rays BA and BC.
See https://en.wikipedia.org/wiki/Angle.

Uses the two-argument form of `atan` returning value in range [-π, π].
See https://en.wikipedia.org/wiki/Atan2.

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

Uses the two-argument form of `atan` returning value in range [-π, π].
See https://en.wikipedia.org/wiki/Atan2.

## Examples

```julia
∠(Vec(1,0), Vec(0,1)) == π/2
```
"""
∠(u::V, v::V) where {V<:Vec{2}} = atan(u × v, u ⋅ v)  # preserve sign
∠(u::V, v::V) where {V<:Vec{3}} = atan(norm(u × v), u ⋅ v)  # discard sign
