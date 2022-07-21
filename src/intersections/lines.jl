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

# compute the parameters of the lines defined by the points a -> b and c -> d 
# so that the resulting points have minimal distance
# the ranks help to identify the different types of intersection
function intersectparameters(a::Point{N,T}, b::Point{N,T}, 
                             c::Point{N,T}, d::Point{N,T}) where {N,T}
  # solves the equation (approximately):
  # a + λ₁ ⋅ v⃗₁ = c + λ₂ ⋅ v⃗₂,
  # with v⃗₁ = b - a, v⃗₂ = d - c
  A = [(b - a) (c - d)]
  y = c - a
  QR = qr([A y])

  # calculate the rank of the augmented matrix
  rₐ = sum(sum(abs.(QR.R), dims = 2) .> atol(T))
  # calculate the rank of the rectangular matrix
  r = sum(abs.(diag(QR.R)[1:2]) .> atol(T))

  if r >= 2 # calculate parameters of intersection or closest point
    QR = qr(A)
    λ = QR.R \ QR.Q' * y
    return λ[1], λ[2], r, rₐ
  else # parallel or collinear
    return zero(T), zero(T), r, rₐ
  end

  # Intersection r == rₐ == 2
  # Collinear: r == rₐ == 1
  # no intersection: r != rₐ
  #   no intersection and parallel:  r == 1, rₐ == 2
  #   no intersection, skew lines: r == 2, rₐ == 3
end
