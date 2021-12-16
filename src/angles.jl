# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ∠(A, B, C)

Angle ∠ABC between rays BA and BC.
See https://en.wikipedia.org/wiki/Angle.

Uses the two-argument form of `atan`.
See https://en.wikipedia.org/wiki/Atan2.

## Examples

```julia
∠(Point(1,0), Point(0,0), Point(0,1)) == π/2
```
"""
function ∠(A::P, B::P, C::P) where {P<:Point{2}}
  BA, BC = A - B, C - B
  atan(BA × BC, BA ⋅ BC) # preserve sign
end

function ∠(A::P, B::P, C::P) where {P<:Point{3}}
  BA, BC = A - B, C - B
  atan(norm(BA × BC), BA ⋅ BC) # discard sign
end

"""
    ∠(u, v)

Angle between vectors u and v.
See https://en.wikipedia.org/wiki/Angle.

Uses [`acos`](@ref) returning a value in the range [0, π].

## Examples

```julia
∠(Vec(1,0), Vec(0,1)) == π/2
```
"""
function ∠(u::Vec{Dim,T}, v::Vec{Dim,T}) where {Dim,T}
  cosθ = (u ⋅ v) / (norm(u) * norm(v))
  acos(clamp(cosθ, -one(T), one(T)))
end
