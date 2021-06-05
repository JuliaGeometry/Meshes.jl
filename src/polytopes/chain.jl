# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Chain(p1, p2, ..., pn)

A polygonal chain from a sequence of points `p1`, `p2`, ..., `pn`.
See https://en.wikipedia.org/wiki/Polygonal_chain.
"""
struct Chain{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Polytope{1,Dim,T}
  vertices::V
end

# specialize Polytope outer constructor to use standard vector
Chain(vertices::Vararg{<:Point}) = Chain(collect(vertices))

Chain(vertices::CircularVector) = Chain([collect(vertices); first(vertices)])

"""
    npoints(chain)

Return the total number of points used to represent the
chain, no matter if it is closed or open.

See also [`nvertices`](@ref).

### Notes

This function is provided for IO purposes. Most algorithms
should be written in terms of `nvertices` and `vertices`
as they are consistent with each other.
"""
npoints(c::Chain) = length(c.vertices)

"""
    nvertices(chain)

Return the number of vertices of the `chain`. In the case
that the chain is closed, the number of vertices is one
less the total number of points used to represent the chain.
"""
nvertices(c::Chain) = npoints(c) - isclosed(c)

"""
    vertices(chain)

Return the vertices of the `chain`. In the case that the
chain is closed, the returned vector is circular and can
be indexed with arbitrarily negative or positive indices.
"""
function vertices(c::Chain)
  if isclosed(c)
    vs = @view c.vertices[begin:end-1]
    CircularVector(vs)
  else
    c.vertices
  end
end

"""
    segments(chain)

Return the segments linking consecutive points of the `chain`.
"""
function segments(c::Chain)
  v = c.vertices
  n = length(v)
  (Segment(view(v, [i,i+1])) for i in 1:n-1)
end

"""
    isclosed(chain)

Tells whether or not the chain is closed.

A closed chain is also known as a ring.
"""
isclosed(c::Chain) = first(c.vertices) == last(c.vertices)

"""
   issimple(chain)

Tells whether or not the `chain` is simple.

A chain is simple when all its segments only
intersect at end points.
"""
function issimple(c::Chain)
  ss = collect(segments(c))
  for i in 1:length(ss)
    for j in i+1:length(ss)
      I = intersecttype(ss[i], ss[j])
      if !(I isa CornerTouchingSegments || I isa NoIntersection)
        return false
      end
    end
  end
  true
end

"""
    windingnumber(point, chain)

Winding number of `point` with respect to the `chain`.
The winding number is the total number of times that
the chain travels counterclockwise around the point.
See https://en.wikipedia.org/wiki/Winding_number.

## References

* Balbes, R. and Siegel, J. 1990. [A robust method for calculating
  the simplicity and orientation of planar polygons]
  (https://www.sciencedirect.com/science/article/abs/pii/0167839691900198)
"""
function windingnumber(p::Point{Dim,T}, c::Chain{Dim,T}) where {Dim,T}
  vₒ, vs = p, c.vertices
  ∑ = sum(∠(vs[i], vₒ, vs[i+1]) for i in 1:length(vs)-1)
  ∑ / T(2π)
end

abstract type OrientationMethod end

struct WindingOrientation <: OrientationMethod end

struct TriangleOrientation <: OrientationMethod end

"""
    orientation(chain, [method])

Returns the orientation of the `chain` as either
counter-clockwise (CCW) or clockwise (CW).

Optionally, specify the orientation `method`:

* `WindingOrientation()` - Balbes, R. and Siegel, J. 1990.
* `TriangleOrientation()` - Held, M. 1998.

Default method is `WindingOrientation()`.

## References

* Balbes, R. and Siegel, J. 1990. [A robust method for calculating
  the simplicity and orientation of planar polygons]
  (https://www.sciencedirect.com/science/article/abs/pii/0167839691900198)
* Held, M. 1998. [FIST: Fast Industrial-Strength Triangulation of Polygons]
  (https://link.springer.com/article/10.1007/s00453-001-0028-4)
"""
orientation(c::Chain) = orientation(c, WindingOrientation())

function orientation(c::Chain{Dim,T}, ::WindingOrientation) where {Dim,T}
  # pick any segment
  x1, x2 = c.vertices[1:2]
  x̄ = centroid(Segment(x1, x2))
  w = T(2π)*windingnumber(x̄, c) - ∠(x1, x̄, x2)
  isapprox(w, T(π), atol=atol(T)) ? :CCW : :CW
end

function orientation(c::Chain{Dim,T}, ::TriangleOrientation) where {Dim,T}
  v = vertices(c)
  Δ(i) = signarea(v[1], v[i], v[i+1])
  a = mapreduce(Δ, +, 2:length(v)-1)
  a ≥ zero(T) ? :CCW : :CW
end

"""
    unique(chain)

Return a new `chain` without duplicate vertices.
"""
Base.unique(c::Chain) = unique!(deepcopy(c))

"""
    unique!(chain)

Remove duplicate vertices in the `chain`.
"""
function Base.unique!(c::Chain)
  # sort vertices lexicographically
  verts = c.vertices
  perms = sortperm(coordinates.(verts))

  # remove true duplicates
  keep = Int[]
  sorted = @view verts[perms]
  for i in 1:length(sorted)-1
    if sorted[i] != sorted[i+1]
      # save index in the original vector
      push!(keep, perms[i])
    end
  end
  push!(keep, last(perms))

  # preserve chain order
  sort!(keep)

  # update vertices in place
  copy!(verts, verts[keep])

  c
end

"""
    close!(chain)

Close the `chain`, i.e. repeat the first vertex
at the end of the vertex list.
"""
function close!(c::Chain)
  push!(c.vertices, first(c.vertices))
  c
end

"""
    open!(chain)

Open the `chain`, i.e. remove the last vertex.
"""
function open!(c::Chain)
  pop!(c.vertices)
  c
end

"""
    reverse!(chain)

Reverse the `chain` vertices in place.
"""
function Base.reverse!(c::Chain)
  reverse!(c.vertices)
  c
end

"""
    reverse(chain)

Reverse the `chain` vertices.
"""
Base.reverse(c::Chain) = reverse!(deepcopy(c))

"""
    angles(chain)

Return angles `∠(vᵢ-₁, vᵢ, vᵢ+₁)` at all vertices
`vᵢ` of the `chain`. If the chain is open, the first
and last vertices have no angles. Positive angles
represent a CCW rotation whereas negative angles
represent a CW rotation. In either case, the
absolute value of the angles returned is never
greater than `π`.
"""
function angles(c::Chain)
  θs = map(2:length(c.vertices)-1) do i
    p1 = c.vertices[i-1]
    p2 = c.vertices[i  ]
    p3 = c.vertices[i+1]
    ∠(p1, p2, p3)
  end

  if isclosed(c)
    p1 = c.vertices[end-1]
    p2 = c.vertices[1]
    p3 = c.vertices[2]
    pushfirst!(θs, ∠(p1, p2, p3))
  end

  θs
end

"""
    innerangles(chain)

Return inner angles of the *closed* `chain`. Inner
angles are always positive, and unlike `angles`
they can be greater than `π`.
"""
function innerangles(c::Chain{Dim,T}) where {Dim,T}
  @assert isclosed(c) "Inner angles only defined for closed chains"

  # correct sign of angles in case orientation is CW
  θs = orientation(c) == :CW ? -angles(c) : angles(c)

  [θ > 0 ? 2*T(π) - θ : -θ for θ in θs]
end

Base.view(c::Chain, inds) = Chain(view(vertices(c), inds))

function Base.show(io::IO, c::Chain{Dim,T}) where {Dim,T}
  N = npoints(c)
  print(io, "$N-Chain{$Dim,$T}")
end

function Base.show(io::IO, ::MIME"text/plain", c::Chain{Dim,T}) where {Dim,T}
  v = c.vertices
  println(io, c)
  print(io, io_lines(v, "  "))
end
