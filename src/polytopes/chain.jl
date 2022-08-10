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
Chain(vertices::Vararg{P}) where {P<:Point} = Chain(collect(vertices))

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
    boundary(chain)

Return the boundary of the `chain`.
"""
function boundary(c::Chain)
  if isclosed(c)
    nothing
  else
    vs = c.vertices
    bs = [first(vs), last(vs)]
    PointSet(bs)
  end
end

"""
    isclosed(chain)

Tells whether or not the chain is closed.

A closed chain is also known as a ring.
"""
isclosed(c::Chain) = first(c.vertices) == last(c.vertices)

"""
    isperiodic(chain)

Tells whether or not the `chain` is periodic
along each parametric dimension.
"""
isperiodic(c::Chain) = (isclosed(c),)

"""
   issimple(chain)

Tells whether or not the `chain` is simple.

A chain is simple when all its segments only
intersect at end points.
"""
function issimple(c::Chain)
  λ(I) = !(type(I) == CornerTouchingSegments ||
           type(I) == NoIntersection)
  ss = collect(segments(c))
  for i in 1:length(ss)
    for j in i+1:length(ss)
      if intersection(λ, ss[i], ss[j])
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
    unique(chain)

Return a new `chain` without duplicate vertices.
"""
Base.unique(c::Chain) = unique!(deepcopy(c))

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
    close(chain)

Same as [`close!`](@ref) but operates on a copy
of the original `chain`.
"""
Base.close(c::Chain) = close!(deepcopy(c))

"""
    open!(chain)

Open the `chain`, i.e. remove the last vertex.
"""
function open!(c::Chain)
  pop!(c.vertices)
  c
end

"""
    open(chain)

Same as [`open!`](@ref) but operates on a copy
of the original `chain`.
"""
Base.open(c::Chain) = open!(deepcopy(c))

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

"""
    bridge(chains; width=0)

Build bridges of given `width` between `chains` of a polygon.

### Notes

- Please read the docstring of the corresponding method for
  [`Polygon`](@ref) for additional details and references.
"""
function bridge(chains::AbstractVector{<:Chain{2,T}}; width=zero(T)) where {T}
  # retrieve chains as vectors of coordinates
  pchains = [coordinates.(vertices(open(c))) for c in chains]

  # sort vertices lexicographically
  coords  = [coord for pchain in pchains for coord in pchain]
  indices = sortperm(sortperm(coords))

  # each chain has its own set of indices
  pinds = Vector{Int}[]; offset = 0
  for nvertex in length.(pchains)
    push!(pinds, indices[offset+1:offset+nvertex])
    offset += nvertex
  end

  # sort chains based on leftmost vertex
  leftmost = argmin.(pinds)
  minimums = getindex.(pinds, leftmost)
  reorder  = sortperm(minimums)
  leftmost = leftmost[reorder]
  minimums = minimums[reorder]
  pchains  = pchains[reorder]
  pinds    = pinds[reorder]

  # initialize outer boundary
  outer = first(pchains)
  oinds = first(pinds)

  # merge holes into outer boundary
  for i in 2:length(pchains)
    inner = pchains[i]
    iinds = pinds[i]
    l = leftmost[i]
    m = minimums[i]

    # find closest vertex in boundary
    dmin, jmin = Inf, 0
    for j in findall(oinds .≤ m)
      d = sum(abs, outer[j] - inner[l])
      if d < dmin
        dmin, jmin = d, j
      end
    end

    # create a bridge of given width δ
    # from line segment A--B. The point
    # A is split into A′ and A′′ and the
    # point B is split into B′ and B′′
    A = outer[jmin]
    B = inner[l]
    δ = width
    v = B - A
    u = Vec(-v[2], v[1])
    n = u / norm(u)
    A′  = A + δ/2 * n
    A′′ = A - δ/2 * n
    B′  = B + δ/2 * n
    B′′ = B - δ/2 * n

    # insert hole at closest vertex
    outer = [
      outer[begin:jmin-1]; [A′, B′];
      circshift(inner, -l+1)[2:end];
      [B′′, A′′]; outer[jmin+1:end]
    ]
    oinds = [
      oinds[begin:jmin];
      circshift(iinds, -l+1);
      [iinds[l]];
      oinds[jmin:end]
    ]
  end

  # find duplicate vertices
  duplicates = Tuple{Int,Int}[]
  occurred = Dict{Int,Int}()
  for (i, ind) in enumerate(oinds)
    if haskey(occurred, ind)
      push!(duplicates, (occurred[ind], i))
    else
      occurred[ind] = i
    end
  end

  # close outer boundary
  push!(outer, first(outer))

  outerchain = Chain(Point.(outer))

  outerchain, duplicates
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
