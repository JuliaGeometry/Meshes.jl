# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# Comparisons between unitful and non-unitful quantities

isequalzero(x) = x == zero(x)
isequalone(x) = x == oneunit(x)

isapproxequal(x, y; kwargs...) = isapprox(x, y; kwargs...)
isapproxzero(x; kwargs...) = isapprox(x, zero(x); kwargs...)
isapproxone(x; kwargs...) = isapprox(x, oneunit(x); kwargs...)

ispositive(x) = x > zero(x)
isnegative(x) = x < zero(x)
isnonpositive(x) = x ≤ zero(x)
isnonnegative(x) = x ≥ zero(x)

"""
    mayberound(λ, x, tol)

Round `λ` to `x` if it is within the tolerance `tol`.
"""
function mayberound(λ, x; kwargs...)
  isapprox(λ, x; kwargs...) ? x : λ
end
