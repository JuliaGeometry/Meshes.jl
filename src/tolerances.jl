# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

const ATOL_FLOAT64 = @load_preference("atol_float64", 1.0e-10)
const ATOL_FLOAT32 = Float32(@load_preference("atol_float32", 1.0f-5))

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
atol(::Type{Float64}) = ATOL_FLOAT64
atol(::Type{Float32}) = ATOL_FLOAT32
atol(ℒ::Type{<:Len}) = atol(numtype(ℒ)) * unit(ℒ)
atol(𝒜::Type{<:Area}) = atol(numtype(𝒜))^2 * unit(𝒜)
atol(𝒱::Type{<:Vol}) = atol(numtype(𝒱))^3 * unit(𝒱)
