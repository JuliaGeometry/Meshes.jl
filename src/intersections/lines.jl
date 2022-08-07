# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

#=
The intersection type can be one of three types:

1. intersect at one point
2. overlap at more than one point
3. do not overlap nor intersect
=#
function intersection(f, l1::Line{N,T}, l2::Line{N,T}) where {N,T}
  a, b = l1(0), l1(1)
  c, d = l2(0), l2(1)

  λ₁, _, r, rₐ = intersectparameters(a, b, c, d)

  if r == rₐ == 2
    return @IT CrossingLines (a + λ₁ * (b - a)) f
  elseif r == rₐ == 1
    return @IT OverlappingLines l1 f
  else
    return @IT NoIntersection nothing f
  end
end

# compute the intersection of two lines assuming that it is a point
function intersectpoint(l1::Line, l2::Line)
  a, b = l1(0), l1(1)
  c, d = l2(0), l2(1)
  λ₁, _ = intersectparameters(a, b, c, d)
  a + λ₁ * (b - a)
 end

 """
    intersectparameters(a, b, c, d)

Compute the parameters `λ₁` and `λ₂` of the lines 
`a + λ₁ ⋅ v⃗₁`, with ` v⃗₁ = b - a` and
`c + λ₂ ⋅ v⃗₂`, with ` v⃗₂ = d - c`
spanned by the input points `a`, `b` resp. `c`, `d`
such that to yield line points with minimal distance 
or the intersection point (if lines intersect).

Furthermore the ranks `r` of the matrix of the linear system 
`A ⋅ λ⃗ = y⃗`, with `A = [v⃗₁ -v⃗₂], y⃗ = c - a`
and the rank `rₐ` of the augmented matrix `[A y⃗]` are calculated in order to 
identify the intersection type:

- Intersection: r == rₐ == 2
- Collinear: r == rₐ == 1
- No intersection: r != rₐ
  - No intersection and parallel:  r == 1, rₐ == 2
  - No intersection, skew lines: r == 2, rₐ == 3
"""
function intersectparameters(a::Point{N,T}, b::Point{N,T}, 
                             c::Point{N,T}, d::Point{N,T}) where {N,T}
  A = [(b - a) (c - d)]
  y = c - a
  QRₐ = qr([A y])

  # calculate the rank of the augmented matrix
  # by checking the zero entries of the diagonal of R
  # for N == 2 one has to check the L1 norm of rows as 
  # there are more columns than rows
  rₐ = N == 2 ? sum(sum(abs, QRₐ.R; dims = 2) .> atol(T)) :
                sum(abs.(diag(QRₐ.R)) .> atol(T))
  # calculate the rank of the rectangular matrix
  r = sum(sum(abs, QRₐ.R[:,1:2]; dims = 2) .> atol(T))

  if r ≥ 2 # calculate parameters of intersection or closest point
    QR = qr(A)
    λ = QR.R \ QR.Q' * y
  else # parallel or collinear
    λ = [zero(T), zero(T)]
  end
  λ[1], λ[2], r, rₐ

end
