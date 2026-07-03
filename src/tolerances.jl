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
atol(::Type{T}) where {T<:AbstractFloat} = eps(T)^(3 // 4)
atol(𝒬::Type{<:Quantity}) = atol(numtype(𝒬)) * unit(𝒬)
atol(𝒟::Type{<:ForwardDiff.Dual}) = atol(ForwardDiff.valtype(𝒟))

# relative tolerance for single and double precision
const RTOL32 = ScopedValue(eps(Float32)^(1 // 3))
const RTOL64 = ScopedValue(eps(Float64)^(1 // 3))

"""
    rtol(T)
    rtol(x::T)

Relative tolerance used in algorithms for approximate
comparison with numbers of type `T`. It is used in the
source code for numerical integration for example.
"""
rtol(x) = rtol(typeof(x))
rtol(::Type{Float32}) = RTOL32[]
rtol(::Type{Float64}) = RTOL64[]
rtol(::Type{T}) where {T<:AbstractFloat} = eps(T)^(1 // 3)
rtol(𝒬::Type{<:Quantity}) = rtol(numtype(𝒬))
rtol(𝒟::Type{<:ForwardDiff.Dual}) = rtol(ForwardDiff.valtype(𝒟))

# maximum length for discretization of non-Euclidean geometries
const MAXLEN = ScopedValue(500u"km")

"""
    maxlen()

Maximum length used for refinement of non-Euclidean geometries.
"""
maxlen() = MAXLEN[]
