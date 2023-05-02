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

isperiodic(::Type{<:Box{Dim}}) where {Dim} = ntuple(i->false, Dim)

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

nvertices(b::Box{Dim}) where Dim = 2^Dim

vertices(b::Box) = collect(vertex(b, ind) for ind in 1:nvertices(b))

function vertex(b::Box{Dim}, ind) where Dim
  1 <= ind <= nvertices(b) ||
    throw(ArgumentError("attempted to access vertex $ind of $(typeof(b))"))
  xmin, xmax = coordinates.(extrema(b))
  ind -= 1 # zero-based index
  coords = ntuple(Dim) do d
    low = iszero(ind & (1 << (d-1)))
    if d == 1 && !iszero(ind & 2)
      low = !low
    end
    low ? xmin[d] : xmax[d]
  end
  Point(coords)
end

function boundary(b::Box{2})
  v = vertices(b)
  Chain([v; first(v)])
end

function boundary(b::Box{3})
  v = vertices(b)
  I = [(4,3,2,1),(6,5,1,2),(3,7,6,2),
       (4,8,7,3),(1,5,8,4),(6,7,8,5)]
  SimpleMesh(v, connect.(I))
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

function (b::Box{Dim,T})(uv...) where {Dim,T}
  if !all(x -> zero(T) ≤ x ≤ one(T), uv)
    throw(DomainError(uv, "b(u, v, ...) is not defined for u, v, ... outside [0, 1]ⁿ."))
  end
  b.min + uv .* (b.max - b.min)
end