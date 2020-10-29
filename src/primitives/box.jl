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
struct Box{N,T} <: Primitive{N,T}
  min::Point{N,T}
  max::Point{N,T}
end

Box(min::Tuple, max::Tuple) = Box(Point(min), Point(max))

Base.minimum(b::Box) = b.min
Base.maximum(b::Box) = b.max
Base.extrema(b::Box) = b.min, b.max
sides(b::Box) = b.max - b.min
volume(b::Box) = prod(b.max - b.min)

function Base.in(p::Point{Dim}, b::Box{Dim}) where {Dim}
  l, u = coordinates.((b.min, b.max))
  x = coordinates(p)
  for i in 1:Dim
    l[i] ≤ x[i] && x[i] ≤ u[i] || return false
  end
  true
end
