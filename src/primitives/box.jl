# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Box(min, max)

An axis-aligned box with `min` and `max` corners.
See https://en.wikipedia.org/wiki/Hyperrectangle.

## Example

```julia
Box(Point(0,0,0), Point(1,1,1)) # unit cube
```
"""
struct Box{Dim,T} <: Primitive{Dim,T}
  min::Point{Dim,T}
  max::Point{Dim,T}
end

Box(min::Tuple, max::Tuple) = Box(Point(min), Point(max))

paramdim(::Type{<:Box{Dim}}) where {Dim} = Dim

isconvex(::Type{<:Box}) = true

Base.minimum(b::Box) = b.min
Base.maximum(b::Box) = b.max
Base.extrema(b::Box) = b.min, b.max
center(b::Box) = Point((coordinates(b.max) + coordinates(b.min)) / 2)
measure(b::Box) = prod(b.max - b.min)
diagonal(b::Box) = norm(b.max - b.min)
sides(b::Box) = b.max - b.min

function vertices(b::Box{2})
  A = coordinates(b.min)
  B = coordinates(b.max)
  Point.([
    (A[1], A[2]),
    (B[1], A[2]),
    (B[1], B[2]),
    (A[1], B[2]),
  ])
end

function vertices(b::Box{3})
  A = coordinates(b.min)
  B = coordinates(b.max)
  Point.([
    (A[1], A[2], A[3]),
    (B[1], A[2], A[3]),
    (B[1], B[2], A[3]),
    (A[1], B[2], A[3]),
    (A[1], A[2], B[3]),
    (B[1], A[2], B[3]),
    (B[1], B[2], B[3]),
    (A[1], B[2], B[3]),
  ])
end

function Base.in(p::Point{Dim}, b::Box{Dim}) where {Dim}
  l, u = coordinates.((b.min, b.max))
  x = coordinates(p)
  for i in 1:Dim
    l[i] ≤ x[i] && x[i] ≤ u[i] || return false
  end
  true
end

Base.issubset(b1::Box{Dim}, b2::Box{Dim}) where {Dim} =
  b1.min ∈ b2 && b1.max ∈ b2
