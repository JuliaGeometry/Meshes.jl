# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Affine(A, b)

Affine transform `Ax + b` with matrix `A` and vector `b`.

# Examples

```julia
Affine(AngleAxis(0.2, 1.0, 0.0, 0.0), [-2, 2, 2])
Affine(Angle2d(π / 2), SVector(2, -2))
Affine([0 -1; 1 0], [-2, 2])
```
"""
struct Affine{Dim,M<:StaticMatrix{Dim,Dim},V<:StaticVector{Dim}} <: CoordinateTransform
  A::M
  b::V
end

function Affine(A::AbstractMatrix, b::AbstractVector)
  sz = size(A)
  if !allequal(sz)
    throw(ArgumentError("`A` must be a square matrix"))
  end
  Dim = first(sz)
  if Dim ≠ length(b)
    throw(ArgumentError("`A` and `b` must have the same dimension"))
  end
  Affine(_assmatrix(Dim, A), _assvector(Dim, b))
end

parameters(t::Affine) = (A=t.A, b=t.b)

isaffine(::Type{<:Affine}) = true

isrevertible(t::Affine) = isinvertible(t)

function isinvertible(t::Affine)
  d = det(t.A)
  !isapprox(d, zero(d), atol=atol(typeof(d)))
end

function inverse(t::Affine)
  A = inv(t.A)
  b = A * -t.b
  Affine(A, b)
end

applycoord(t::Affine, v::Vec) = t.A * v

applycoord(t::Affine, p::Point) = Point(t.A * coordinates(p) + t.b)

# --------------
# SPECIAL CASES
# --------------

applycoord(t::Affine, b::Box{2}) = applycoord(t, convert(Quadrangle, b))

applycoord(t::Affine, b::Box{3}) = applycoord(t, convert(Hexahedron, b))

# -----------------
# HELPER FUNCTIONS
# -----------------

_assmatrix(Dim, A::StaticMatrix) = A
_assmatrix(Dim, A) = SMatrix{Dim,Dim}(A)

_assvector(Dim, b::StaticVector) = b
_assvector(Dim, b) = SVector{Dim}(b)
