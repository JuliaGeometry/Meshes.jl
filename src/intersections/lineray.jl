# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

#=
The intersection type can be one of four types:

1. intersect at one inner pont (CrossingLineRay -> Point)
2. intersect at origin of ray (MidTouchingLineRay -> Point)
3. overlap at more than one point (OverlappingLineRay -> Ray)
5. do not overlap nor intersect (NoIntersection)
=#
function intersection(f::Function, l1::Line{3,T}, r1::Ray{3,T}) where {T}
  a, b = l1(0), l1(1)
  c, d = r1(0), r1(1)

  if measure(Tetrahedron(a, b, c, d)) > 0 # not in same plane
    return @IT NoIntersection nothing f
  elseif isapprox(norm((b - a) × (c - d)), zero(T), atol=atol(T)^2) #parallel 
    if isapprox(measure(Triangle(a, b, c)), zero(T), atol=atol(T)^2) # collinear
        return @IT OverlappingLineRay (r1) f #OverlappingRaySegment
	else # parallel with distance > 0
    	return @IT NoIntersection nothing f
	end
  else # vectors in same plane but not parallel
    v1  = a - b
  	v2  = c - d
  	v12 = a - c 
	s = (v12[1] * v2[2] - v12[2] * v2[1]) / (v1[1] * v2[2] - v1[2] * v2[1])
	if s < 0
		return @IT NoIntersection nothing f
	else
        if r1.p ∈ l1
            return @IT MidTouchingLineRay (r1.p) f
        else    
		    return @IT CrossingLineRay (a-s*v1) f # returns point
        end
	end
    
  end
end
#=
The intersection type can be one of four types:

1. intersect at one inner pont (CrossingLineRay -> Point)
2. intersect at origin of ray (MidTouchingLineRay -> Point)
3. overlap at more than one point (OverlappingLineRay -> Ray)
5. do not overlap nor intersect (NoIntersection)
=#
function intersection(f::Function, l1::Line{2,T}, r1::Ray{2,T}) where {T}
    a, b = l1(0), l1(1)
    c, d = r1(0), r1(1)
  
   if isapprox(norm((b - a) × (c - d)), zero(T), atol=atol(T)^2) #parallel 
      if isapprox(measure(Triangle(a, b, c)), zero(T), atol=atol(T)^2) # collinear
          return @IT OverlappingLineRay (r1) f #OverlappingRaySegment
      else # parallel with distance > 0
          return @IT NoIntersection nothing f
      end
    else # vectors in same plane but not parallel
      v1  = a - b
        v2  = c - d
        v12 = a - c 
      s = (v12[1] * v2[2] - v12[2] * v2[1]) / (v1[1] * v2[2] - v1[2] * v2[1])
      if s < 0
          return @IT NoIntersection nothing f
      else
          if r1.p ∈ l1
              return @IT MidTouchingLineRay (r1.p) f
          else    
              return @IT CrossingLineRay (a-s*v1) f # returns point
          end
      end
      
    end
  end

# for 2D and 3D use lines.jl implementation
# NOTE: no check whether resulting point is in ray
function intersectpoint(r1::Ray, s1::Line)
   return inersectpoint(Line(r1(0), r1(1)), Line(s1(0), s1(1)))
end
