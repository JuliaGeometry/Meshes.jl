# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ray(p, v)

A ray originating at point `p`, pointed in direction `v`.
It can be called as `r(t)` with `t > 0` to cast it at
`p + t * v`.
"""
struct Ray{Dim,T} <: Primitive{Dim,T}
  p::Point{Dim,T}
  v::Vec{Dim,T}
end

Ray(p::Tuple, v::Tuple) = Ray(Point(p), Vec(v))

paramdim(::Type{<:Ray}) = 1

isconvex(::Type{<:Ray}) = true

boundary(r::Ray) = r.p

function (r::Ray)(t)
  if t < 0
    throw(DomainError(t, "r(t) is not defined for t < 0."))
  end
  r.p + t * r.v
end

"""
    origin(ray)

The starting point of the ray.
"""
origin(r::Ray) = r.p

"""
    direction(ray)

The direction of the ray.
"""
direction(r::Ray) = r.v


function Base.in(p::Point, r::Ray)
	w = norm(r.v)
	d = evaluate(Euclidean(), p, Line(r(0), r(1)))
	if d+w ≈ w
		return sum((p - r.p) .* r.v) >= 0 # return true if point is in ray
	else
		return false
	end
end

==(r1::Ray, r2::Ray) = r1.p ∈ r2 && r2.p ∈ r1 && r1(1) ∈ r2 && r2(1) ∈ r1 
