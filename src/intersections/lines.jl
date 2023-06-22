# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

#=
The intersection type can be one of three types:

1. intersect at one point
2. overlap at more than one point
3. do not overlap nor intersect
=#
function intersection(f, line₁::Line, line₂::Line)
  a, b = line₁(0), line₁(1)
  c, d = line₂(0), line₂(1)

  λ₁, _, r, rₐ = intersectparameters(a, b, c, d)

  if r == rₐ == 2
    return @IT CrossingLines (a + λ₁ * (b - a)) f
  elseif r == rₐ == 1
    return @IT OverlappingLines line₁ f
  else
    return @IT NoIntersection nothing f
  end
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
- Colinear: r == rₐ == 1
- No intersection: r != rₐ
  - No intersection and parallel:  r == 1, rₐ == 2
  - No intersection, skew lines: r == 2, rₐ == 3
"""
function intersectparameters(a::Point{Dim,T}, b::Point{Dim,T}, c::Point{Dim,T}, d::Point{Dim,T}) where {Dim,T}
  A = [(b - a) (c - d)]
  y = c - a

  # calculate the rank of the augmented matrix by checking
  # the zero entries of the diagonal of R
  _, R = qr([A y])

  # for Dim == 2 one has to check the L1 norm of rows as 
  # there are more columns than rows
  rₐ = sum(sum(abs, R, dims=2) .> atol(T))

  # calculate the rank of the rectangular matrix
  r = sum(sum(abs, R[:, 1:2], dims=2) .> atol(T))

  # calculate parameters of intersection or closest point
  if r ≥ 2
    λ = A \ y
  else # parallel or collinear
    λ = SVector(zero(T), zero(T))
  end

  λ[1], λ[2], r, rₐ
end
