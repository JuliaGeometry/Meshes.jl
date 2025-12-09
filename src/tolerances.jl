# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# absolute tolerance for single and double precision
const ATOL32 = ScopedValue(eps(Float32)^(3 // 4))
const ATOL64 = ScopedValue(eps(Float64)^(3 // 4))

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
atol(::Type{Float32}) = ATOL32[]
atol(::Type{Float64}) = ATOL64[]
atol(ℒ::Type{<:Len}) = atol(numtype(ℒ)) * unit(ℒ)
atol(𝒜::Type{<:Area}) = atol(numtype(𝒜))^2 * unit(𝒜)
atol(𝒱::Type{<:Vol}) = atol(numtype(𝒱))^3 * unit(𝒱)

# maximum length for discretization of non-Euclidean geometries
const MAXLEN = ScopedValue(500u"km")

"""
Maximum length used for discretization of non-Euclidean geometries.
It is used in the source code in calls to the [`discretize`](@ref)
and [`simplexify`](@ref) functions.
"""
maxlen() = MAXLEN[]
