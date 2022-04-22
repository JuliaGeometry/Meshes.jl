# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

#=
The intersection type can be one of three types:

1. intersect at one point
2. overlap at more than one point
3. do not overlap nor intersect
=#
function intersection(f::Function, l1::Line{3,T}, l2::Line{3,T}) where {T}
  a, b = l1(0), l1(1)
  c, d = l2(0), l2(1)

  if measure(Tetrahedron(a, b, c, d)) > 0 # not in same plane
    return @IT NoIntersection nothing f
  elseif isapprox(norm((b - a) × (c - d)), zero(T), atol=atol(T)^2)
    if a in l2
      return @IT OverlappingLines l1 f
    else # parallel lines
      return @IT NoIntersection nothing f
    end
  else
    return @IT CrossingLines intersectpoint(l1, l2) f
  end
end

function intersection(f::Function, l1::Line{2,T}, l2::Line{2,T}) where {T}
  a, b = l1(0), l1(1)
  c, d = l2(0), l2(1)

  if !isapprox(abs((b - a) × (c - d)), zero(T), atol=atol(T)^2)
    return @IT CrossingLines intersectpoint(l1, l2) f
  elseif isapprox(measure(Triangle(a, b, c)), zero(T), atol=atol(T)^2)
    return @IT OverlappingLines l1 f
  else
    return @IT NoIntersection nothing f
  end
end

# compute the intersection of two lines assuming that it is a point
function intersectpoint(l1::Line{2}, l2::Line{2})
  a, b = l1(0), l1(1)
  c, d = l2(0), l2(1)

  v1  = a - b
  v2  = c - d
  v12 = a - c

  # the intersection point lies in between a and b at a fraction s
  # (https://en.wikipedia.org/wiki/Line-line_intersection#Formulas)
  s = (v12[1] * v2[2] - v12[2] * v2[1]) / (v1[1] * v2[2] - v1[2] * v2[1])

  a - s*v1
end

# compute the intersection of two lines assuming that it is a point 
# same code as for 2D case should also work in 3D by ignoring third component
# TODO: find a better method
function intersectpoint(l1::Line{3}, l2::Line{3})
  a, b = l1(0), l1(1)
  c, d = l2(0), l2(1)

  v1  = a - b
  v2  = c - d
  v12 = a - c

  # the intersection point lies in between a and b at a fraction s
  # (https://en.wikipedia.org/wiki/Line-line_intersection#Formulas)
  s = (v12[1] * v2[2] - v12[2] * v2[1]) / (v1[1] * v2[2] - v1[2] * v2[1])

  a - s*v1
end