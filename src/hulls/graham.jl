# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    GrahamScan

Compute the convex hull of a set of points or geometries using the
Graham's scan method. See [https://en.wikipedia.org/wiki/Graham_scan]
(https://en.wikipedia.org/wiki/Graham_scan).

The method has complexity `O(n*log(n))` where `n` is the number of points.

## References

* Cormen et al. 2009. [Introduction to Algorithms]
  (https://mitpress.mit.edu/books/introduction-algorithms-third-edition)
"""
struct GrahamScan <: HullMethod end

function hull(pset::PointSet{2,T}, ::GrahamScan) where {T}
  # remove duplicates
  Q = coordinates.(pset) |> unique

  # sort by y then by x
  sort!(Q, by=reverse)

  # corner cases
  n = length(Q)
  n == 1 && return Point(Q[1])
  n == 2 && return Segment(Point(Q[1]), Point(Q[2]))
  if n == 3
    p₀, p₁, p₂ = Q
    θ = ∠(Point(p₁), Point(p₀), Point(p₂))
    if isapprox(θ, zero(T), atol=atol(T))
      return Segment(Point(p₀), Point(p₂))
    else
      c = Chain(Point(p₀), Point(p₁), Point(p₂), Point(p₀))
      c = orientation(c) == :CCW ? c : reverse(c)
      return PolyArea(c)
    end
  end

  # sort by polar angle
  p₀ = Point(Q[1])
  p  = Point.(Q[2:n])
  x  = p₀ + Vec{2,T}(1, 0)
  θ  = [∠(x, p₀, pᵢ) for pᵢ in p]
  p  = p[sortperm(θ)]

  # rotational sweep
  c = [p₀, p[1], p[2]]
  for pᵢ in p[3:end]
    while ∠(c[end-1], c[end], pᵢ) > atol(T)
      pop!(c)
    end
    push!(c, pᵢ)
  end

  # close chain
  push!(c, c[1])

  PolyArea(c)
end


"""
  hull_simple(p, h, s, v)

calculate the convex hull for the matrix of points.  This is actually an implementation of the pseudocode from  [https://en.wikipedia.org/wiki/Graham_scan]  (https://en.wikipedia.org/wiki/Graham_scan).

This is a fast and convinient method when point coordinates are in a matrix, and you have to create 10000000 hulls in a second.  

parameters:   
 * p is the matrix [y x w]  of the points; where 'y' is the point dimension and w is the number of points.  
 * h is the preallocated vector of the indices of the points in p [w]
 * s is the preallocated vector of the indices of the points in p, which makes a hull   [w]
 * v preallocated vector of Folat64 [w]

return the number of the points in a hull, and indices of the hull points in 's'
"""
function hull_simple(p, h, s, v)
	j3, n = size(p)
	p0 = 1

	# detect p0
	@inbounds for i=2:n
		if p[2, i] < p[2, p0]   # minimum y
			p0 = i
		elseif p[2, i] == p[2, p0]
			if p[1, i] < p[1, p0]
				p0 = i
			end
		end
	end
	#@info "p0 = $p0"
	
	@inbounds for i = 1:n
		if i == p0
			v[i] = 1.0
			continue
		end
		v[i] = (((p[1, i] - p[1, p0])) / sqrt((p[1, i] - p[1, p0])^2 + (p[2, i] - p[2, p0])^2)) # this is the cosinus of the angle
	end
	#@info "not sorted v = $v"
	sortperm!(h, v, rev = true)			
	#@info "sorted v = $v"
	u = 0
	@inbounds for i = 1:n
		# v = (p2[1] - p1[1])*(p3[2] - p1[2]) - (p2[2] - p1[2])*(p3[1] - p1[1])
		while u > 1 &&  (p[1, s[u]] - p[1, s[u-1]])*(p[2, h[i]] - p[2, s[u-1]]) - (p[2, s[u]] - p[2, s[u-1]])*(p[1, h[i]] - p[1, s[u-1]]) <= 0.0
			u -= 1
		end
		u += 1
		s[u] = h[i]
	end
	return u
end


"""
	zeroInside(p, s)

check if zero point [0.0, 0.0] is inside the convex hull (hull created by hull_simple() function). 

 * p is the matrix [y x w]  of the points; where 'y' is the point dimension and w is the number of points.  
 * s is the vector of indices in p, of the points that we use for convex hull
 * u number of the indices in 's'

 return true if zero point is inside the hull.
"""
function zeroInside(p, s, u)
	for i = 2:u
		if (p[1, s[i]] - p[1, s[i-1]]) * (-p[2, s[i-1]]) - (p[2, s[i]] - p[2, s[i-1]])*(-p[1, s[i-1]]) < 0.0
			return false
		end
	end
	# the last segment:
	if (p[1, s[1]] - p[1, s[u]])*(-p[2, s[u]]) - (p[2, s[1]] - p[2, s[u]])*(-p[1, s[u]]) < 0.0
		return false
	end
	return true
end

export hull_simple, zeroInside