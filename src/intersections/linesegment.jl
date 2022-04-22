# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

#=
The intersection type can be one of four types:

1. intersect at one inner point (CrossingLineSegment -> Point)
2. intersect at one corner point of segment (MidTouchingLineSegment -> Point)
3. overlap at more than one point (OverlappingLineSegment -> Segment)
4. do not overlap nor intersect (NoIntersection)
=#
function intersection(f::Function, l1::Line{3,T}, s1::Segment{3,T}) where {T}
  a, b = l1(0), l1(1)
  c, d = s1(0), s1(1)

  if measure(Tetrahedron(a, b, c, d)) > 0 # not in same plane
    return @IT NoIntersection nothing f
  elseif isapprox(norm((b - a) × (c - d)), zero(T), atol=atol(T)^2) #parallel 
    if isapprox(measure(Triangle(a, b, c)), zero(T), atol=atol(T)^2) # collinear
        return @IT OverlappingLineSegment s1 f
	else # parallel with distance > 0
    	return @IT NoIntersection nothing f
	end
  else # vectors in same plane but not parallel
    v1  = a - b
    v2  = c - d
    v12 = a - c 
    #TODO: find better way to calculate parameters in 3D including third component
    s = (v12[1] * v1[2] - v12[2] * v1[1]) / (v1[1] * v2[2] - v1[2] * v2[1])
    if s < 0 || s > 1 
        return @IT NoIntersection nothing f
    elseif isapprox(s, zero(T), atol=atol(T)) || isapprox(s, one(T), atol=atol(T)) 
        return @IT MidTouchingLineSegment (s < 0.5 ? c : d) f #MidTouchingLineSegment
    else
        return @IT CrossingLineSegment (c-s*v2) f #CrossingRaySegment # returns point # equal to c - s * v2
    end
    
  end
end

#=
The intersection type can be one of four types:

1. intersect at one inner point (CrossingLineSegment -> Point)
2. intersect at one corner point of segment (MidTouchingLineSegment -> Point)
3. overlap at more than one point (OverlappingLineSegment -> Segment)
4. do not overlap nor intersect (NoIntersection)
=#
function intersection(f::Function, l1::Line{2,T}, s1::Segment{2,T}) where {T}
    a, b = l1(0), l1(1)
    c, d = s1(0), s1(1)
  
    if isapprox(norm((b - a) × (c - d)), zero(T), atol=atol(T)^2) #parallel 
      if isapprox(measure(Triangle(a, b, c)), zero(T), atol=atol(T)^2) # collinear
          return @IT OverlappingLineSegment s1 f
      else # parallel with distance > 0
          return @IT NoIntersection nothing f
      end
    else # vectors in same plane but not parallel
      v1  = a - b
      v2  = c - d
      v12 = a - c 
      #TODO: find better way to calculate parameters in 3D including third component
      s = (v12[1] * v1[2] - v12[2] * v1[1]) / (v1[1] * v2[2] - v1[2] * v2[1])
      if s < 0 || s > 1 
          return @IT NoIntersection nothing f
      elseif isapprox(s, zero(T), atol=atol(T)) || isapprox(s, one(T), atol=atol(T)) 
          return @IT MidTouchingLineSegment (s < 0.5 ? c : d) f #MidTouchingLineSegment
      else
          return @IT CrossingLineSegment (c-s*v2) f #CrossingRaySegment # returns point # equal to c - s * v2
      end
      
    end
  end

# for 2D and 3D use lines.jl implementation
# NOTE: no check whether resulting point is in ray and segment
function intersectpoint(l1::Line, s1::Segment)
   return inersectpoint(Line(l1(0), l1(1)), Line(s1(0), s1(1)))
end
