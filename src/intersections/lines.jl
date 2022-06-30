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

  p1, _, r, rₐ = intersectparameters(a,b,c,d)

  if r == rₐ == 2
    return @IT CrossingLines a + p1 * (b-a) f
  elseif r == rₐ == 1
    return @IT OverlappingLines l1 f
  else
    @IT NoIntersection nothing f
  end

  #= 
  if !isapprox(measure(Tetrahedron(a, b, c, d)), zero(T), atol = atol(T)) # not in same plane
    return @IT NoIntersection nothing f
  elseif isapprox(norm((b - a) × (c - d)), zero(T), atol=atol(T))
    if a in l2 # collinear
      return @IT OverlappingLines l1 f
    else # parallel lines
      return @IT NoIntersection nothing f
    end
  else
    return @IT CrossingLines intersectpoint(l1, l2) f
  end
  =#
end

# compute the intersection of two lines assuming that it is a point
function intersectpoint(l1::Line, l2::Line)
  a, b = l1(0), l1(1)
  c, d = l2(0), l2(1)
  p1, _ = intersectparameters(a,b,c,d)
  a + p1 * (b - a)
 end

# compute the intersection parameters of the lines defined by the points a -> b and c -> d 
function intersectparameters_fast(a::Point, b::Point, c::Point, d::Point) 
  v1  = a - b
  v2  = c - d
  v12 = a - c
  
  # the intersection point lies in between a and b at a fraction p1
  # respectively between c and d at a fraction p2
  # (https://en.wikipedia.org/wiki/Line-line_intersection#Formulas)
  p1 = (v12[1] * v2[2] - v12[2] * v2[1]) / (v1[1] * v2[2] - v1[2] * v2[1])
  p2 = (v12[1] * v1[2] - v12[2] * v1[1]) / (v1[1] * v2[2] - v1[2] * v2[1])
  p1, p2
end

# compute the parameters of the lines defined by the points a -> b and c -> d 
# so that the resulting points have minimal distance
# the ranks help to identify the different types of intersection

function intersectparameters(a::Point{N,T}, b::Point{N,T}, c::Point{N,T}, d::Point{N,T}) where {N,T}
  # solves the equation (approximately): a + λ₁ ⋅ v⃗₁ = c + λ₂ ⋅ v⃗₂, v⃗₁ = b - a, v⃗₂ = d - c
  A = [(b - a) (c - d)]
  b = c - a

  λ = A \ b

  # calculate the rank of the augmented matrix
  rₐ = rank([A b], atol = atol(T))
  # calculate the rank of the rectangular matrix
  r = rank(A, atol = atol(T))

  # Intersection r == rₐ == 2
  # Collinear: r == rₐ == 1
  # no intersection: r != rₐ
    # no intersection and parallel:  r == 1, rₐ == 2
    # no intersection, skew lines: r == 2, rₐ == 3
	λ[1], λ[2], r, rₐ
end
