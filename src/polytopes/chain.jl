# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Chain(p1, p2, ..., pn)

A polygonal chain from a sequence of points `p1`, `p2`, ..., `pn`.
See https://en.wikipedia.org/wiki/Polygonal_chain.

See also [`Rope`](@ref) and [`Ring`](@ref).
"""
abstract type Chain{Dim,T} <: Polytope{1,Dim,T} end

"""
    Rope(p1, p2, ..., pn)

An open polygonal chain from a sequence of points `p1`, `p2`, ..., `pn`.

See also [`Chain`](@ref) and [`Ring`](@ref).
"""
struct Rope{Dim,T,V<:AbstractVector{Point{Dim,T}}} <: Chain{Dim,T}
  vertices::V
end

"""
    Ring(p1, p2, ..., pn)

A closed polygonal chain from a sequence of points `p1`, `p2`, ..., `pn`.

See also [`Chain`](@ref) and [`Rope`](@ref).
"""
struct Ring{Dim,T,V<:CircularVector{Point{Dim,T}}} <: Chain{Dim,T}
  vertices::V
end

Ring(vertices::AbstractVector{P}) where {P<:Point} = Ring(CircularVector(vertices))

"""
    segments(chain)

Return the segments linking consecutive points of the `chain`.
"""
function segments(c::Chain)
  v = c.vertices
  n = length(v) - !isclosed(c)
  (Segment(view(v, [i, i + 1])) for i in 1:n)
end

"""
    boundary(chain)

Return the boundary of the `chain`.
"""
function boundary(r::Rope)
  v = r.vertices
  PointSet([first(v), last(v)])
end
boundary(::Ring) = nothing

"""
    isclosed(chain)

Tells whether or not the chain is closed.

A closed chain is also known as a ring.
"""
isclosed(c::Rope) = false
isclosed(c::Ring) = true

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
  λ(I) = !(type(I) == CornerTouchingSegments || type(I) == NoIntersection)
  ss = collect(segments(c))
  for i in 1:length(ss)
    for j in (i + 1):length(ss)
      if intersection(λ, ss[i], ss[j])
        return false
      end
    end
  end
  true
end

"""
    windingnumber(point, ring)

Winding number of `point` with respect to the `ring`.
The winding number is the total number of times that
the ring travels counterclockwise around the point.
See https://en.wikipedia.org/wiki/Winding_number.

## References

* Balbes, R. and Siegel, J. 1990. [A robust method for calculating
  the simplicity and orientation of planar polygons]
  (https://www.sciencedirect.com/science/article/abs/pii/0167839691900198)
"""
function windingnumber(p::Point{2,T}, r::Ring{2,T}) where {T}
  v = r.vertices
  n = length(v)
  ∑ = sum(∠(v[i], p, v[i + 1]) for i in 1:n)
  ∑ / T(2π)
end

abstract type OrientationMethod end

struct WindingOrientation <: OrientationMethod end

struct TriangleOrientation <: OrientationMethod end

"""
    orientation(ring, [method])

Returns the orientation of the `ring` as either
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
orientation(r::Ring) = orientation(r, WindingOrientation())

function orientation(r::Ring{2,T}, ::WindingOrientation) where {T}
  # pick any segment
  x1, x2 = r.vertices[1:2]
  x̄ = centroid(Segment(x1, x2))
  w = T(2π) * windingnumber(x̄, r) - ∠(x1, x̄, x2)
  isapprox(w, T(π), atol=atol(T)) ? :CCW : :CW
end

function orientation(r::Ring{2,T}, ::TriangleOrientation) where {T}
  v = vertices(r)
  Δ(i) = signarea(v[1], v[i], v[i + 1])
  a = mapreduce(Δ, +, 2:(length(v) - 1))
  a ≥ zero(T) ? :CCW : :CW
end

"""
    unique!(chain)

Remove duplicate vertices in the `chain`.
Closed chains remain closed.
"""
function Base.unique!(c::Chain)
  # sort vertices lexicographically
  verts = open(c).vertices # work with underlying array, even for closed chains
  perms = sortperm(coordinates.(verts))

  # remove true duplicates
  keep = Int[]
  sorted = @view verts[perms]
  for i in 1:(length(sorted) - 1)
    if sorted[i] != sorted[i + 1]
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
Closed chains remain closed.
"""
Base.unique(c::Chain) = unique!(deepcopy(c))

"""
    close(chain)

Close the `chain`, i.e. add a segment going from the last to the first vertex.
"""
Base.close(r::Rope) = Ring(r.vertices)
Base.close(r::Ring) = r

"""
    open(chain)

Open the `chain`, i.e. remove the segment going from the last to the first vertex.
"""
# call `open` again to avoid issues in case of nested CircularVector
Base.open(r::Rope) = r
Base.open(r::Ring) = open(Rope(parent(r.vertices)))

"""
    reverse!(chain)

Reverse the `chain` vertices in place.
"""
# do not change which vertex comes first for closed chains
Base.reverse!(r::Rope) = (reverse!(r.vertices); r)
Base.reverse!(r::Ring) = (reverse!(@view r.vertices[(begin + 1):end]); r)

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
  vs = vertices(c)
  i1 = firstindex(vs) + !isclosed(c)
  i2 = lastindex(vs) - !isclosed(c)
  map(i -> ∠(vs[i - 1], vs[i], vs[i + 1]), i1:i2)
end

"""
    innerangles(ring)

Return inner angles of the `ring`. Inner
angles are always positive, and unlike
`angles` they can be greater than `π`.
"""
function innerangles(r::Ring{2,T}) where {T}
  # correct sign of angles in case orientation is CW
  θs = orientation(r) == :CW ? -angles(r) : angles(r)
  [θ > 0 ? 2 * T(π) - θ : -θ for θ in θs]
end

"""
    bridge(rings; width=0)

Build bridges of given `width` between `rings` of a polygon.

### Notes

- Please read the docstring of the corresponding method for
  [`Polygon`](@ref) for additional details and references.
"""
function bridge(rings::AbstractVector{<:Ring{2,T}}; width=zero(T)) where {T}
  # retrieve chains as vectors of coordinates
  pchains = [coordinates.(vertices(open(r))) for r in rings]

  # sort vertices lexicographically
  coords = [coord for pchain in pchains for coord in pchain]
  indices = sortperm(sortperm(coords))

  # each chain has its own set of indices
  pinds = Vector{Int}[]
  offset = 0
  for nvertex in length.(pchains)
    push!(pinds, indices[(offset + 1):(offset + nvertex)])
    offset += nvertex
  end

  # sort chains based on leftmost vertex
  leftmost = argmin.(pinds)
  minimums = getindex.(pinds, leftmost)
  reorder = sortperm(minimums)
  leftmost = leftmost[reorder]
  minimums = minimums[reorder]
  pchains = pchains[reorder]
  pinds = pinds[reorder]

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
    dmin, jmin = typemax(T), 0
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
    A′ = A + δ / 2 * n
    A′′ = A - δ / 2 * n
    B′ = B + δ / 2 * n
    B′′ = B - δ / 2 * n

    # insert hole at closest vertex
    outer = [
      outer[begin:(jmin - 1)]
      [A′, B′]
      circshift(inner, -l + 1)[2:end]
      [B′′, A′′]
      outer[(jmin + 1):end]
    ]
    oinds = [
      oinds[begin:jmin]
      circshift(iinds, -l + 1)
      [iinds[l]]
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
  outerring = Ring(Point.(outer))

  outerring, duplicates
end

function Base.show(io::IO, c::Chain{Dim,T}) where {Dim,T}
  n = nvertices(c)
  name = nameof(typeof(c))
  print(io, "$n-$name{$Dim,$T}")
end

function Base.show(io::IO, ::MIME"text/plain", c::Chain{Dim,T}) where {Dim,T}
  println(io, c)
  print(io, io_lines(c.vertices, "  "))
end
