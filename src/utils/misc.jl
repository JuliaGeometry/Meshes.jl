# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    signarea(A, B, C)

Compute signed area of triangle formed by points `A`, `B` and `C`.
"""
function signarea(A::Point, B::Point, C::Point)
  assertion(embeddim(A) == 2, "points must be 2-dimensional")
  ((B - A) × (C - A)) / 2
end

"""
    householderbasis(n)

Returns a pair of orthonormal tangent vectors `u` and `v` from a normal `n`,
such that `u`, `v`, and `n` form a right-hand orthogonal system.

## References

* D.S. Lopes et al. 2013. ["Tangent vectors to a 3-D surface normal: A geometric tool
  to find orthogonal vectors based on the Householder transformation"]
  (https://doi.org/10.1016/j.cad.2012.11.003)
"""
function householderbasis(n::Vec{3,ℒ}) where {ℒ}
  n̂ = norm(n)
  i = argmax(n .+ n̂)
  n̂ᵢ = Vec(ntuple(j -> j == i ? n̂ : zero(ℒ), 3))
  h = n + n̂ᵢ
  H = (I - 2h * transpose(h) / (transpose(h) * h)) * unit(ℒ)
  u, v = [H[:, j] for j in 1:3 if j != i]
  i == 2 && ((u, v) = (v, u))
  Vec(u), Vec(v)
end

"""
    svdbasis(points)

Returns the 2D basis that retains most of the variance in the list of 3D `points`
using the singular value decomposition (SVD).

See <https://math.stackexchange.com/a/99317>.
"""
function svdbasis(p::AbstractVector{<:Point})
  assertion(embeddim(first(p)) == 3, "points must be 3-dimensional")
  ℒ = lentype(eltype(p))
  X = stack(to, p)
  μ = sum(X, dims=2) / size(X, 2)
  Z = X .- μ
  U = usvd(Z).U
  u = Vec(U[:, 1]...)
  v = Vec(U[:, 2]...)
  n = Vec(zero(ℒ), zero(ℒ), oneunit(ℒ))
  isnegative((u × v) ⋅ n) ? (v, u) : (u, v)
end
