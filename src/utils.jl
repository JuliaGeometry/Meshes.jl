# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    signarea(A, B, C)

Compute signed area of triangle formed
by points `A`, `B` and `C`.
"""
function signarea(A::Point{2}, B::Point{2}, C::Point{2})
  ((B - A) × (C - A)) / 2
end

"""
    iscollinear(A, B, C)

Tells whether or not the points
`A`, `B` and `C` are collinear.
"""
function iscollinear(A::Point{Dim,T}, B::Point{Dim,T}, C::Point{Dim,T}) where {Dim,T}
  # points A, B, C are collinear if and only if the
  # cross-products for segments AB and AC with respect
  # to all possible pairs of coordinates are zero
  AB, AC = B - A, C - A
  result = true
  for i in 1:Dim, j in (i+1):Dim
    u = Vec{2,T}(AB[i], AB[j])
    v = Vec{2,T}(AC[i], AC[j])
    if !isapprox(u × v, zero(T), atol=atol(T)^2)
      result = false
      break
    end
  end
  result
end

"""
    sideof(point, segment)

Determines on which side of the oriented `segment`
the `point` lies. Possible results are `:LEFT`,
`:RIGHT` or `:ON` the segment.
"""
function sideof(p::Point{2,T}, s::Segment{2,T}) where {T}
  a, b = vertices(s)
  area = signarea(p, a, b)
  ifelse(area > atol(T), :LEFT, ifelse(area < -atol(T), :RIGHT, :ON))
end

"""
    sideof(point, chain)

Determines on which side of the closed `chain` the
`point` lies. Possible results are `:INSIDE` or
`:OUTSIDE` the chain.
"""
function sideof(p::Point{2,T}, c::Chain{2,T}) where {T}
  w = windingnumber(p, c)
  ifelse(isapprox(w, zero(T), atol=atol(T)), :OUTSIDE, :INSIDE)
end

"""
    mahalanobis(radii, angles, convention)

Return the Mahalanobis distance corresponding to an ellipsoid
with `radii` rotated by given `angles` according to given `convention`.

- For 2D ellipses, there are two radii and one rotation angle.
- For 3D ellipsoids, there are three radii and three rotation angles.

The list of available conventions can be found with:

```julia
julia> subtypes(RotationConvention)
```
"""
function mahalanobis(radii, angles, convention)
  ndims, nangles = length(radii), length(angles)
  valid = (ndims == 3 && nangles == 3) || (ndims == 2 && nangles == 1)
  @assert valid "invalid number of radii/angles"

  # invert radii if necessary
  invert = mainaxis(convention) == :Y
  ranges = invert ? [radii[i] for i in reverse(1:ndims,1,2)] : vec(radii)

  # scaling matrix
  λ = one(eltype(ranges)) ./ ranges.^2
  Λ = Diagonal(SVector{ndims}(λ))

  # rotation matrix
  R = rotmat(angles, convention)

  # ellipsoid metric
  W = Symmetric(R'*Λ*R)
  Mahalanobis(W)
end