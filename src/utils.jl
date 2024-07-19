# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# auxiliary type for dispatch purposes
const GeometryOrDomain = Union{Geometry,Domain}

"""
    assertion(cond, msg)

Throws an `AssertionError(msg)` if `cond` is `false`.
"""
assertion(cond, msg) = cond || throw(AssertionError(msg))

"""
    assertdim(geom, dim)

Throws an `AssertionError` if the `embeddim` of the geometry `geom`
is different than the specified dimension `dim`. 
"""
assertdim(geom, dim) = assertion(embeddim(geom) == dim, lazy"geometry must be embedded in $dim-dimensional space")

"""
    constructor(G)

Given a (parametric) type `G{T₁,T₂,...}`, return the type `G`.
"""
constructor(G::Type{<:GeometryOrDomain}) = getfield(Meshes, nameof(G))

"""
    fitdims(dims, D)

Fit tuple `dims` to a given length `D` by repeating the last dimension.
"""
function fitdims(dims::Dims{N}, D) where {N}
  ntuple(i -> i ≤ N ? dims[i] : last(dims), D)
end

"""
    collectat(iter, inds)

Collect iterator `iter` at indices `inds` without materialization.
"""
function collectat(iter, inds)
  if isempty(inds)
    eltype(iter)[]
  else
    selectat(inds) = enumerate ⨟ TakeWhile(x -> first(x) ≤ last(inds)) ⨟ Filter(y -> first(y) ∈ inds) ⨟ Map(last)
    iter |> selectat(inds) |> tcollect
  end
end

"""
    withcrs(g, v)

Point at the end of the vector `v` with the same CRS of `g`.
"""
function withcrs(g::GeometryOrDomain, v::StaticVector)
  C = crs(g)
  cart = Cartesian{datum(C)}(Tuple(v))
  ctor = CoordRefSystems.constructor(C)
  Point(convert(ctor, cart))
end

"""
    flat(p)

Flatten coordinates of point `p` to Cartesian coordinates,
ignoring the original units of the coordinate reference system.
"""
flat(p::Point) = Point(flat(coords(p)))
flat(c::LatLon) = Cartesian{datum(c)}(CoordRefSystems.rawvalues(c))
flat(c::CRS) = convert(Cartesian, c)

"""
    signarea(A, B, C)

Compute signed area of triangle formed by points `A`, `B` and `C`.
"""
function signarea(A::Point, B::Point, C::Point)
  assertdim(A, 2)
  ((B - A) × (C - A)) / 2
end

"""
    householderbasis(n)

Returns a pair of orthonormal tangent vectors `u` and `v` from a normal `n`,
such that `u`, `v`, and `n` form a right-hand orthogonal system.

## References

* D.S. Lopes et al. 2013. ["Tangent vectors to a 3-D surface normal: A geometric tool
  to find orthogonal vectors based on the Householder transformation"]
  (https://doi.org/10.1016/j.cad.2012.11.003)
"""
function householderbasis(n::Vec{3,ℒ}) where {ℒ}
  n̂ = norm(n)
  i = argmax(n .+ n̂)
  n̂ᵢ = Vec(ntuple(j -> j == i ? n̂ : zero(ℒ), 3))
  h = n + n̂ᵢ
  H = (I - 2h * transpose(h) / (transpose(h) * h)) * unit(ℒ)
  u, v = [H[:, j] for j in 1:3 if j != i]
  i == 2 && ((u, v) = (v, u))
  Vec(u), Vec(v)
end

"""
    svdbasis(points)

Returns the 2D basis that retains most of the variance in the list of 3D `points`
using the singular value decomposition (SVD).

See <https://math.stackexchange.com/a/99317>.
"""
function svdbasis(p::AbstractVector{<:Point})
  assertdim(first(p), 3)
  ℒ = lentype(eltype(p))
  X = reduce(hcat, to.(p))
  μ = sum(X, dims=2) / size(X, 2)
  Z = X .- μ
  U = usvd(Z).U
  u = Vec(U[:, 1]...)
  v = Vec(U[:, 2]...)
  n = Vec(zero(ℒ), zero(ℒ), oneunit(ℒ))
  isnegative((u × v) ⋅ n) ? (v, u) : (u, v)
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
function intersectparameters(a::Point, b::Point, c::Point, d::Point)
  A = ustrip.([(b - a) (c - d)])
  y = ustrip.(c - a)
  T = eltype(A)

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

"""
    XYZ(xyz)

Generate the coordinate arrays `XYZ` from the coordinate vectors `xyz`.
"""
@generated function XYZ(xyz::NTuple{Dim,<:AbstractVector{T}}) where {Dim,T}
  exprs = ntuple(Dim) do d
    quote
      a = xyz[$d]
      A = Array{T,Dim}(undef, length.(xyz))
      @nloops $Dim i A begin
        @nref($Dim, A, i) = a[$(Symbol(:i_, d))]
      end
      A
    end
  end
  Expr(:tuple, exprs...)
end

isequalzero(x) = x == zero(x)
isequalone(x) = x == oneunit(x)

isapproxequal(x, y; atol=atol(x), kwargs...) = isapprox(x, y; atol, kwargs...)
isapproxzero(x; atol=atol(x), kwargs...) = isapprox(x, zero(x); atol, kwargs...)
isapproxone(x; atol=atol(x), kwargs...) = isapprox(x, oneunit(x); atol, kwargs...)

ispositive(x) = x > zero(x)
isnegative(x) = x < zero(x)
isnonpositive(x) = x ≤ zero(x)
isnonnegative(x) = x ≥ zero(x)

# Function wrappers that handle units
# The result units of some operations, such as dot and cross, 
# are treated in a special way to handle Meshes.jl use cases

function usvd(A)
  u = unit(eltype(A))
  F = svd(ustrip.(A))
  SVD(F.U * u, F.S * u, F.Vt * u)
end

uinv(A) = inv(ustrip.(A)) * unit(eltype(A))^-1

unormalize(a::Vec{Dim,ℒ}) where {Dim,ℒ} = Vec(normalize(a) * unit(ℒ))

udot(a::Vec{Dim,ℒ}, b::Vec{Dim,ℒ}) where {Dim,ℒ} = ustrip(a ⋅ b) * unit(ℒ)
udot(a::Vec{Dim,ℒ₁}, b::Vec{Dim,ℒ₂}) where {Dim,ℒ₁,ℒ₂} = udot(promote(a, b)...)

ucross(a::Vec{Dim,ℒ}, b::Vec{Dim,ℒ}) where {Dim,ℒ} = Vec(ustrip.(a × b) * unit(ℒ))
ucross(a::Vec{Dim,ℒ₁}, b::Vec{Dim,ℒ₂}) where {Dim,ℒ₁,ℒ₂} = ucross(promote(a, b)...)

ucross(a::Vec{Dim,ℒ}, b::Vec{Dim,ℒ}, c::Vec{Dim,ℒ}) where {Dim,ℒ} = Vec(ustrip.(a × b × c) * unit(ℒ))

urotbetween(u::Vec, v::Vec) = rotation_between(ustrip.(u), ustrip.(v))

urotapply(R::Rotation, v::Vec{Dim,ℒ}) where {Dim,ℒ} = Vec(R * ustrip.(v) * unit(ℒ))
