# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const ATOL64 = ScopedValue(1e-10)
const ATOL32 = ScopedValue(1f-5)

"""
    atol(T)
    atol(x::T)

Absolute tolerance used in algorithms for approximate
comparison with numbers of type `T`. It is used in the
source code in calls to the [`isapprox`](@ref) function:

```julia
isapprox(a::T, b::T, atol=atol(T))
```
"""
atol(x) = atol(typeof(x))
@inline atol(::Type{Float64}) = ATOL64[]::Float64
@inline atol(::Type{Float32}) = ATOL32[]::Float32
atol(ℒ::Type{<:Len}) = atol(numtype(ℒ)) * unit(ℒ)
atol(𝒜::Type{<:Area}) = atol(numtype(𝒜))^2 * unit(𝒜)
atol(𝒱::Type{<:Vol}) = atol(numtype(𝒱))^3 * unit(𝒱)
