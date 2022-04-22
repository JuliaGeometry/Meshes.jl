# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

#=
The intersection type can be one of five types:

1. intersect at one inner point (CrossingRaySegment -> Point)
2. intersect at one corner point of segment xor origin of ray (MidTouchingRaySegment -> Point)
3. intersects at one corner point of segment and origin of ray (CornerTouchingRaySegment -> Point)
4. overlap at more than one point (OverlappingRaySegment -> Segment)
5. do not overlap nor intersect (NoIntersection)
=#
function intersection(f::Function, r1::Ray{3,T}, s1::Segment{3,T}) where {T}
  a, b = r1(0), r1(1)
  c, d = s1(0), s1(1)

  if measure(Tetrahedron(a, b, c, d)) > 0 # not in same plane
    return @IT NoIntersection nothing f
  elseif isapprox(norm((b - a) × (c - d)), zero(T), atol=atol(T)^2) #parallel 
    if isapprox(measure(Triangle(a, b, c)), zero(T), atol=atol(T)^2) # collinear
        r_c = mean((c - a) ./ r1.v) # has to be the same
		r_d = mean((d - a) ./ r1.v)
		if r_c > 0 
			if r_d > 0
				return @IT OverlappingRaySegment s1 f# OverlappingRaySegment
			else
				return @IT OverlappingRaySegment Segment(r1.p, c) f# OverlappingRaySegment
			end
		else
			if r_d > 0
				return @IT OverlappingRaySegment (Segment(r1.p, d)) f #OverlappingRaySegment
			else
				return @IT NoIntersection nothing f
			end
		end
	else # parallel with distance > 0
    	return @IT NoIntersection nothing f
	end
  else # vectors in same plane but not parallel
    v1  = a - b
    v2  = c - d
    v12 = a - c 
    #TODO: find better way to calculate parameters in 3D including third component
    r = (v12[1] * v2[2] - v12[2] * v2[1]) / (v1[1] * v2[2] - v1[2] * v2[1])
    s = (v12[1] * v1[2] - v12[2] * v1[1]) / (v1[1] * v2[2] - v1[2] * v2[1])
    if r < 0 || (s < 0 || s > 1)
        
        return @IT NoIntersection nothing f
    elseif isapprox(r, zero(T), atol=atol(T))
        if isapprox(s, one(T), atol=atol(T)) || isapprox(s,zero(T), atol=atol(T))
            return @IT CornerTouchingRaySegment a f #CornerTouchingRaySegment
        else
            return @IT MidTouchingRaySegment a f #MidTouchingRaySegment
        end
    else
        if isapprox(s, one(T), atol=atol(T)) || isapprox(s,zero(T), atol=atol(T))
            return @IT MidTouchingRaySegment (s < 0.5 ? c : d)  f# MidTouchingRaySegment
        else
            return @IT CrossingRaySegment (a-r*v1) f #CrossingRaySegment # returns point # equal to c - s * v2
        end
    end
    
  end
end

#=
The intersection type can be one of five types:

1. intersect at one inner point (CrossingRaySegment -> Point)
2. intersect at one corner point of segment xor origin of ray (MidTouchingRaySegment -> Point)
3. intersects at one corner point of segment and origin of ray (CornerTouchingRaySegment -> Point)
4. overlap at more than one point (OverlappingRaySegment -> Segment)
5. do not overlap nor intersect (NoIntersection)
=#
function Meshes.intersection(f::Function, r1::Ray{2,T}, s1::Segment{2,T}) where {T}
	a, b = r1(0), r1(1)
	c, d = s1(0), s1(1)

	if !isapprox(abs((b - a) × (c - d)), zero(T), atol=atol(T)^2) # not parallel
		
		v1  = a - b
  		v2  = c - d
  		v12 = a - c 
		r = (v12[1] * v2[2] - v12[2] * v2[1]) / (v1[1] * v2[2] - v1[2] * v2[1])
		s = (v12[1] * v1[2] - v12[2] * v1[1]) / (v1[1] * v2[2] - v1[2] * v2[1])
		if r < 0 || (s < 0 || s > 1)
			
			return @IT NoIntersection nothing f
		elseif isapprox(r, zero(T), atol=atol(T))
			if isapprox(s, one(T), atol=atol(T)) || isapprox(s,zero(T), atol=atol(T))
				return @IT CornerTouchingRaySegment a f #CornerTouchingRaySegment
			else
				return @IT MidTouchingRaySegment a f #MidTouchingRaySegment
			end
		else
			if isapprox(s, one(T), atol=atol(T)) || isapprox(s,zero(T), atol=atol(T))
				return @IT MidTouchingRaySegment (s < 0.5 ? c : d)  f# MidTouchingRaySegment
			else
				return @IT CrossingRaySegment (a-r*v1) f #CrossingRaySegment # returns point # equal to c - s * v2
			end
		end
	elseif isapprox(measure(Triangle(a, b, c)), zero(T), atol=atol(T)^2) # parallel and even overlapping
		# identify segment points c and d via ray parameter r_c and r_d 
        # (at least one has to be positive)
		r_c = mean((c - a) ./ r1.v) # has to be the same
		r_d = mean((d - a) ./ r1.v)
		if r_c > 0 
			if r_d > 0
				return @IT OverlappingRaySegment s1 f# OverlappingRaySegment
			else
				return @IT OverlappingRaySegment Segment(r1.p, c) f# OverlappingRaySegment
			end
		else
			if r_d > 0
				return @IT OverlappingRaySegment (Segment(r1.p, d)) f #OverlappingRaySegment
			else
				return @IT NoIntersection nothing f
			end
		end
	else # parallel with distance > 0
    	return @IT NoIntersection nothing f
	end
end

# for 2D and 3D use lines.jl implementation
# NOTE: no check whether resulting point is in ray and segment
function intersectpoint(r1::Ray, s1::Segment)
   return inersectpoint(Line(r1(0), r1(1)), Line(s1(0), s1(1)))
end
