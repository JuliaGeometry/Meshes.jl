# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    atol(T)

Absolute tolerance used in algorithms for approximate
comparison with numbers of type `T`. It is used in the
source code in calls to the [`isapprox`](@ref) function:

```julia
isapprox(a::T, b::T, atol=atol(T))
```
"""
atol(::Type{Float64}) = 1e-10
atol(::Type{Float32}) = 1f-5
