# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    signarea(A, B, C)

Compute signed area of triangle formed by points `A`, `B` and `C`.
"""
function signarea(A::Point, B::Point, C::Point)
  checkdim(A, 2)
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
  checkdim(first(p), 3)
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

"""
    intersectparameters(a, b, c, d)

Compute the parameters `λ₁` and `λ₂` of the lines 
`a + λ₁ ⋅ v⃗₁`, with `v⃗₁ = b - a` and
`c + λ₂ ⋅ v⃗₂`, with `v⃗₂ = d - c` spanned by the input
points `a`, `b` resp. `c`, `d` such that to yield line
points with minimal distance or the intersection point
(if lines intersect).

Furthermore, the ranks `r` of the matrix of the linear
system `A ⋅ λ⃗ = y⃗`, with `A = [v⃗₁ -v⃗₂], y⃗ = c - a`
and the rank `rₐ` of the augmented matrix `[A y⃗]` are
calculated in order to identify the intersection type:

- Intersection: r == rₐ == 2
- Collinear: r == rₐ == 1
- No intersection: r != rₐ
  - No intersection and parallel:  r == 1, rₐ == 2
  - No intersection, skew lines: r == 2, rₐ == 3
"""
function intersectparameters(a::Point, b::Point, c::Point, d::Point)
  # augmented linear system
  A = ustrip.([(b - a) (c - d) (c - a)])

  # normalize by maximum absolute coordinate
  A = A / maximum(abs, A)

  # numerical tolerance
  T = eltype(A)
  τ = atol(T)

  # check if a vector is non zero
  isnonzero(v) = !isapprox(v, zero(v), atol=τ)

  # calculate ranks by checking the zero rows of
  # the factor R in the QR matrix factorization
  _, R = qr(A)
  r = sum(isnonzero, eachrow(R[:, SVector(1, 2)]))
  rₐ = sum(isnonzero, eachrow(R))

  # calculate parameters of intersection
  if r ≥ 2
    λ = A[:, SVector(1, 2)] \ A[:, 3]
    λ₁, λ₂ = λ[1], λ[2]
  else # parallel or collinear
    λ₁, λ₂ = zero(T), zero(T)
  end

  λ₁, λ₂, r, rₐ
end
