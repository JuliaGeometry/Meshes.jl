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
sides(b::Box) = b.max - b.min
volume(b::Box) = prod(b.max - b.min)

# set operations

"""
    union(b1::Box, b2::Box)

Union between boxes.
"""
function Base.union(b1::Box{N,T1}, b2::Box{N,T2}) where {N,T1,T2}
  m1, M1 = minimum(b1), maximum(b1)
  m2, M2 = minimum(b2), maximum(b2)
  T = promote_type(T1, T2)
  m = min.(coordinates(m1), coordinates(m2))
  M = max.(coordinates(M1), coordinates(M2))
  Box{N,T}(Point(m), Point(M))
end

"""
    intersect(b1::Box, b2::Box)

Intersection between boxes.
"""
function intersect(b1::Box{N,T1}, b2::Box{N,T2}) where {N,T1,T2}
  m1, M1 = minimum(b1), maximum(b1)
  m2, M2 = minimum(b2), maximum(b2)
  T = promote_type(T1, T2)
  m = max.(coordinates(m1), coordinates(m2))
  M = min.(coordinates(M1), coordinates(M2))
  Box{N,T}(Point(m), Point(M))
end

function before(b1::Box{N}, b2::Box{N}) where {N}
  M1 = coordinates(maximum(b1))
  m2 = coordinates(minimum(b2))
  for i in 1:N
    M1[i] < m2[i] || return false
  end
  true
end

function overlaps(b1::Box{N}, b2::Box{N}) where {N}
  m1 = coordinates(minimum(b1))
  M1 = coordinates(maximum(b1))
  m2 = coordinates(minimum(b2))
  M2 = coordinates(maximum(b2))
  for i in 1:N
    M2[i] > M1[i] > m2[i] && m1[i] < m2[i] || return false
  end
  true
end

function starts(b1::Box{N}, b2::Box{N}) where {N}
  m1 = coordinates(minimum(b1))
  M1 = coordinates(maximum(b1))
  m2 = coordinates(minimum(b2))
  M2 = coordinates(maximum(b2))
  if m1 == m2
    for i in 1:N
      M1[i] < M2[i] || return false
    end
    true
  else
    false
  end
end

function during(b1::Box{N}, b2::Box{N}) where {N}
  m1 = coordinates(minimum(b1))
  M1 = coordinates(maximum(b1))
  m2 = coordinates(minimum(b2))
  M2 = coordinates(maximum(b2))
  for i in 1:N
    M1[i] < M2[i] && m1[i] > m2[i] || return false
  end
  true
end

function finishes(b1::Box{N}, b2::Box{N}) where {N}
  m1 = coordinates(minimum(b1))
  M1 = coordinates(maximum(b1))
  m2 = coordinates(minimum(b2))
  M2 = coordinates(maximum(b2))
  if M1 == M2
    for i in 1:N
      m1[i] > m2[i] || return false
    end
    true
  else
    false
  end
end

# containment

"""
    in(b1::Box, b2::Box)

Check if Box `b1` is contained in `b2`. This does not use
strict inequality, so Rects may share faces and this will still
return true.
"""
function Base.in(b1::Box{N}, b2::Box{N}) where {N}
  m1 = coordinates(minimum(b1))
  M1 = coordinates(maximum(b1))
  m2 = coordinates(minimum(b2))
  M2 = coordinates(maximum(b2))
  for i in 1:N
    M1[i] ≤ M2[i] && m1[i] ≥ m2[i] || return false
  end
  true
end

"""
    in(p::Point, b::Box)

Check if a point is in the box.
"""
function Base.in(p::Point{N,T}, b::Box{N,T}) where {T,N}
  m = coordinates(minimum(b))
  M = coordinates(maximum(b))
  x = coordinates(p)
  for i in 1:N
    m[i] ≤ x[i] && x[i] ≤ M[i] || return false
  end
  true
end
