# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# auxiliary type for dispatch purposes
const GeometryOrDomain = Union{Geometry,Domain}

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
    sideof(point, line)

Determines on which side the `point` is in relation to the `line`.
Possible results are `:LEFT`, `:RIGHT` or `:ON` the line.

### Notes

* For a line `l` this method assumes the orientation
  of the line segment formed by the points `l(0)` and `l(1)`.
"""
function sideof(p::Point{2,T}, l::Line{2,T}) where {T}
  area = signarea(p, l(0), l(1))
  ifelse(area > atol(T), :LEFT, ifelse(area < -atol(T), :RIGHT, :ON))
end

"""
    sideof(point, ring)

Determines on which side of the `ring` the `point` lies.
Possible results are `:INSIDE` or `:OUTSIDE` the ring.
"""
function sideof(p::Point{2,T}, r::Ring{2,T}) where {T}
  w = windingnumber(p, r)
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
  r = Ray(point, Vec(zero(T), zero(T), 2 * (z[2] - z[1])))

  hasintersect = false
  edgecrosses = 0
  ps = Point{3,T}[]
  for t in mesh
    I = intersection(r, t)
    if type(I) == Crossing
      hasintersect = !hasintersect
    elseif type(I) ∈ (EdgeTouching, CornerTouching, Touching)
      return :ON
    elseif type(I) == EdgeCrossing
      edgecrosses += 1
    elseif type(I) == CornerCrossing
      p = get(I)
      if !any(≈(p), ps)
        push!(ps, p)
        hasintersect = !hasintersect
      end
    end
  end

  # check how many edges we crossed
  isodd(edgecrosses ÷ 2) && (hasintersect = !hasintersect)
  hasintersect ? (return :INSIDE) : (return :OUTSIDE)
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
  _, i = findmax(n .+ n̂)
  ei = 1:3 .== i
  h = n + n̂ * ei
  H = I - 2h * transpose(h) / (transpose(h) * h)
  u, v = [H[:, j] for j in 1:3 if j != i]
  if i == 2
    u, v = v, u
  end
  u, v
end

"""
    mayberound(λ, x, tol)

Round `λ` to `x` if it is within the tolerance `tol`.
"""
function mayberound(λ::T, x::T, atol=atol(T)) where {T}
  isapprox(λ, x, atol=atol) ? x : λ
end

"""
    intersectparameters(a, b, c, d)

Compute the parameters `λ₁` and `λ₂` of the lines 
`a + λ₁ ⋅ v⃗₁`, with `v⃗₁ = b - a` and
`c + λ₂ ⋅ v⃗₂`, with `v⃗₂ = d - c` spanned by the input
points `a`, `b` resp. `c`, `d` such that to yield line
points with minimal distance or the intersection point
(if lines intersect).

Furthermore, the ranks `r` of the matrix of the linear
system `A ⋅ λ⃗ = y⃗`, with `A = [v⃗₁ -v⃗₂], y⃗ = c - a`
and the rank `rₐ` of the augmented matrix `[A y⃗]` are
calculated in order to identify the intersection type:

- Intersection: r == rₐ == 2
- Colinear: r == rₐ == 1
- No intersection: r != rₐ
  - No intersection and parallel:  r == 1, rₐ == 2
  - No intersection, skew lines: r == 2, rₐ == 3
"""
function intersectparameters(a::Point{Dim,T}, b::Point{Dim,T}, c::Point{Dim,T}, d::Point{Dim,T}) where {Dim,T}
  A = [(b - a) (c - d)]
  y = c - a

  # calculate the rank of the augmented matrix by checking
  # the zero entries of the diagonal of R
  _, R = qr([A y])

  # for Dim == 2 one has to check the L1 norm of rows as 
  # there are more columns than rows
  τ = atol(T)
  rₐ = sum(>(τ), sum(abs, R, dims=2))

  # calculate the rank of the rectangular matrix
  r = sum(>(τ), sum(abs, view(R, :, 1:2), dims=2))

  # calculate parameters of intersection or closest point
  if r ≥ 2
    λ = A \ y
    λ₁, λ₂ = λ[1], λ[2]
  else # parallel or collinear
    λ₁, λ₂ = zero(T), zero(T)
  end

  λ₁, λ₂, r, rₐ
end
