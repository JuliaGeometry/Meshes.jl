# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ring(p1, p2, ..., pn)

A closed polygonal chain from a sequence of points `p1`, `p2`, ..., `pn`.

See also [`Chain`](@ref) and [`Rope`](@ref).
"""
struct Ring{Dim,T,V<:CircularVector{Point{Dim,T}}} <: Chain{Dim,T}
  vertices::V

  function Ring{Dim,T,V}(vertices) where {Dim,T,V}
    if first(vertices) == last(vertices) && length(vertices) ≥ 2
      throw(ArgumentError("""
      First and last vertices of `Ring` constructor must be different
      in the latest version of Meshes.jl. The type itself now holds
      this connectivity information.
      """))
    end
    new(vertices)
  end
end

Ring(vertices::CircularVector{P}) where {P<:Point} = Ring{embeddim(P),coordtype(P),typeof(vertices)}(vertices)
Ring(vertices::AbstractVector{P}) where {P<:Point} = Ring(CircularVector(vertices))

boundary(::Ring) = nothing

isclosed(::Type{<:Ring}) = true

Base.close(r::Ring) = r

# call `open` again to avoid issues in case of nested CircularVector
Base.open(r::Ring) = open(Rope(parent(r.vertices)))

# do not change which vertex comes first for closed chains
Base.reverse!(r::Ring) = (reverse!(@view r.vertices[(begin + 1):end]); r)

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

orientation(r::Ring{3}, method::OrientationMethod) = orientation(Ring(proj2D(vertices(r))), method)

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
