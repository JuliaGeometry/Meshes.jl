# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    fitdims(dims, D)

Fit tuple `dims` to a given length `D` by repeating the last dimension.
"""
function fitdims(dims::Dims{N}, D) where {N}
  ntuple(i -> i ≤ N ? dims[i] : last(dims), D)
end

"""
    signarea(A, B, C)

Compute signed area of triangle formed by points `A`, `B` and `C`.
"""
function signarea(A::Point{2}, B::Point{2}, C::Point{2})
  ((B - A) × (C - A)) / 2
end

"""
    iscollinear(A, B, C)

Tells whether or not the points `A`, `B` and `C` are collinear.
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
    iscoplanar(A, B, C, D)

Tells whether or not the points `A`, `B`, `C` and `D` are coplanar.
"""
function iscoplanar(A::Point{3,T}, B::Point{3,T}, C::Point{3,T}, D::Point{3,T}) where {T}
  vol = volume(Tetrahedron(A, B, C, D))
  isapprox(vol, zero(T), atol = atol(T))
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
    sideof(point, mesh)

Determines whether a `point` is inside, outside or on the surface of a `mesh`.
Possible results are `:INSIDE`, `:OUTSIDE`, or `:ON`.

### Notes

Uses a ray-casting algorithm.
"""
function sideof(point::Point{3,T}, mesh::Mesh{3,T}) where {T}
  @assert paramdim(mesh) == 2 "sideof only defined for surface meshes"
  (eltype(mesh) <: Triangle) || return sideof(point, simplexify(mesh))

  z = last.(coordinates.(extrema(mesh)))
  r = Ray(point, Vec{3,T}(zero(T), zero(T), 2 * (z[2] - z[1])))

  intersects = false
  edgecrosses = 0
  ps = Point{3,T}[]
  for t in mesh
    I = intersection(r, t)
    if type(I) == CrossingRayTriangle
      intersects = !intersects
    elseif type(I) ∈ (EdgeTouchingRayTriangle,
                      CornerTouchingRayTriangle,
                      TouchingRayTriangle)
      return :ON
    elseif type(I) == EdgeCrossingRayTriangle
      edgecrosses += 1
    elseif type(I) == CornerCrossingRayTriangle
      p = get(I)
      if !any(≈(p), ps)
        push!(ps, p)
        intersects = !intersects
      end
    end
  end

  # check how many edges we crossed
  isodd(edgecrosses ÷ 2) && (intersects = !intersects)
  intersects ? (return :INSIDE) : (return :OUTSIDE)
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

Returns the unitless type of a (unitful) type or value. See Unitful.jl.
i.e. `dropunits(1u"mm") == Int`
"""
dropunits(v) = typeof(one(v))

"""
    householderbasis(n)

Returns a pair of orthonormal tangent vectors `u` and `v` from a normal `n`,
such that `u`, `v`, and `n` form a right-hand orthogonal system.

## References

* D.S. Lopes et al. 2013. ["Tangent vectors to a 3-D surface normal: A geometric tool
  to find orthogonal vectors based on the Householder transformation"]
  (https://doi.org/10.1016/j.cad.2012.11.003)
"""
function householderbasis(n)
  n̂ = norm(n)
  _, i = findmax(n.+n̂)
  ei = 1:3 .== i
  h = n + n̂*ei
  H = I - 2h*transpose(h)/(transpose(h)*h)
  u, v  = [H[:,j] for j = 1:3 if j != i]
  if i == 2
      u, v = v, u
  end
  u, v
end

"""
    mayberound(λ, x, tol)

Round `λ` to `x` if it is within the tolerance `tol`.
"""
function mayberound(λ::T, x, atol = atol(T)) where {T}
  isapprox(λ, x, atol = atol) ? x : λ
end

"""
    uvrotation(u, v)

Return rotation matrix from vector `u` to vector `v`.
"""
function uvrotation(u::Vec{2}, v::Vec{2})
  convert(DCM, CounterClockwiseAngle(∠(u, v)))
end

function uvrotation(u::Vec{3}, v::Vec{3})
  u⃗ = normalize(u)
  v⃗ = normalize(v)
  re = √((1 + u⃗ ⋅ v⃗) / 2)
  im = (u⃗ × v⃗) / 2re
  convert(DCM, Quaternion(re, im...))
end
