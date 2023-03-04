# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Box(min, max)

An axis-aligned box with `min` and `max` corners.
See https://en.wikipedia.org/wiki/Hyperrectangle.

## Examples

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

measure(b::Box) = prod(b.max - b.min)

Base.length(b::Box{1}) = measure(b)

area(b::Box{2}) = measure(b)

volume(b::Box{3}) = measure(b)

center(b::Box) = Point((coordinates(b.max) + coordinates(b.min)) / 2)

diagonal(b::Box) = norm(b.max - b.min)

sides(b::Box) = Tuple(b.max - b.min)

vertices(b::Box{1}) = [b.min, b.max]

function vertices(b::Box{2})
  A = coordinates(b.min)
  B = coordinates(b.max)
  return Point.([(A[1], A[2]), (B[1], A[2]), (B[1], B[2]), (A[1], B[2])])
end

function vertices(b::Box{3})
  A = coordinates(b.min)
  B = coordinates(b.max)
  return Point.([
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

function boundary(b::Box{2})
  v = vertices(b)
  return Chain([v; first(v)])
end

function boundary(b::Box{3})
  v = vertices(b)
  I = [(4, 3, 2, 1), (6, 5, 1, 2), (3, 7, 6, 2), (4, 8, 7, 3), (1, 5, 8, 4), (6, 7, 8, 5)]
  return SimpleMesh(v, connect.(I))
end

function Base.in(p::Point{Dim}, b::Box{Dim}) where {Dim}
  l, u = coordinates.((b.min, b.max))
  x = coordinates(p)
  for i in 1:Dim
    l[i] ≤ x[i] && x[i] ≤ u[i] || return false
  end
  return true
end

Base.issubset(b1::Box{Dim}, b2::Box{Dim}) where {Dim} = b1.min ∈ b2 && b1.max ∈ b2
