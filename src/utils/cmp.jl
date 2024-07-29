# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# Comparisons between unitful and non-unitful quantities

isequalzero(x) = x == zero(x)
isequalone(x) = x == oneunit(x)

isapproxequal(x, y; atol=atol(x), kwargs...) = isapprox(x, y; atol, kwargs...)
isapproxzero(x; atol=atol(x), kwargs...) = isapprox(x, zero(x); atol, kwargs...)
isapproxone(x; atol=atol(x), kwargs...) = isapprox(x, oneunit(x); atol, kwargs...)

ispositive(x) = x > zero(x)
isnegative(x) = x < zero(x)
isnonpositive(x) = x ≤ zero(x)
isnonnegative(x) = x ≥ zero(x)

"""
    mayberound(λ, x, tol)

Round `λ` to `x` if it is within the tolerance `tol`.
"""
function mayberound(λ::T, x::T, atol=atol(T)) where {T}
  isapprox(λ, x, atol=atol) ? x : λ
end
