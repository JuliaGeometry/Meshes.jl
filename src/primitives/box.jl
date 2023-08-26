# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Box(min, max)

An axis-aligned box with `min` and `max` corners.
See https://en.wikipedia.org/wiki/Hyperrectangle.

## Examples

```julia
Box(Point(0, 0, 0), Point(1, 1, 1))
Box((0, 0), (1, 1))
```
"""
struct Box{Dim,T} <: Primitive{Dim,T}
  min::Point{Dim,T}
  max::Point{Dim,T}

  function Box{Dim,T}(min, max) where {Dim,T}
    @assert min ⪯ max "`min` must be less than or equal to `max`"
    new(min, max)
  end
end

Box(min::Point{Dim,T}, max::Point{Dim,T}) where {Dim,T} = Box{Dim,T}(min, max)

Box(min::Tuple, max::Tuple) = Box(Point(min), Point(max))

paramdim(::Type{<:Box{Dim}}) where {Dim} = Dim

Base.minimum(b::Box) = b.min

Base.maximum(b::Box) = b.max

Base.extrema(b::Box) = b.min, b.max

center(b::Box) = Point((coordinates(b.max) + coordinates(b.min)) / 2)

diagonal(b::Box) = norm(b.max - b.min)

sides(b::Box) = Tuple(b.max - b.min)

boundary(b::Box{1}) = PointSet([b.min, b.max])

function boundary(b::Box{2})
  A = coordinates(b.min)
  B = coordinates(b.max)
  v = Point.([(A[1], A[2]), (B[1], A[2]), (B[1], B[2]), (A[1], B[2])])
  Ring(v)
end

function boundary(b::Box{3})
  A = coordinates(b.min)
  B = coordinates(b.max)
  v =
    Point.([
      (A[1], A[2], A[3]),
      (B[1], A[2], A[3]),
      (B[1], B[2], A[3]),
      (A[1], B[2], A[3]),
      (A[1], A[2], B[3]),
      (B[1], A[2], B[3]),
      (B[1], B[2], B[3]),
      (A[1], B[2], B[3])
    ])
  c = [(4, 3, 2, 1), (6, 5, 1, 2), (3, 7, 6, 2), (4, 8, 7, 3), (1, 5, 8, 4), (6, 7, 8, 5)]
  SimpleMesh(v, connect.(c))
end

Base.isapprox(b₁::Box, b₂::Box) = b₁.min ≈ b₂.min && b₁.max ≈ b₂.max

function (b::Box{Dim,T})(uv...) where {Dim,T}
  if !all(x -> zero(T) ≤ x ≤ one(T), uv)
    throw(DomainError(uv, "b(u, v, ...) is not defined for u, v, ... outside [0, 1]ⁿ."))
  end
  b.min + uv .* (b.max - b.min)
end

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Box{Dim,T}}) where {Dim,T}
  min = rand(rng, Point{Dim,T})
  max = min + rand(rng, Vec{Dim,T})
  Box(min, max)
end
