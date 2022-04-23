# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    fitdims(dims, D)

Fit tuple `dims` to a given length `D`
by repeating the last dimension.
"""
function fitdims(dims::Dims{N}, D) where {N}
  ntuple(i -> i ≤ N ? dims[i] : last(dims), D)
end

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
    proj2D(points)

Convert a vector of 3D `points` into a vector of 2D
points living in a plane of maximum variance using SVD.
"""
function proj2D(points::AbstractVector{Point{3,T}}) where {T}
  # https://math.stackexchange.com/a/99317
  X = mapreduce(coordinates, hcat, points)
  μ = sum(X, dims=2) / size(X, 2)
  Z = X .- μ
  U = svd(Z).U
  u = U[:,1]
  v = U[:,2]
  [Point(z⋅u, z⋅v) for z in eachcol(Z)]
end

"""
    dropunits(T)

Gets the dimensionless memory type representation of a variable or type.
i.e. `dropunits(1u"mm") == Int`
"""
dropunits(T) = typeof(one(T))
