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

Ring(vertices::CircularVector{Point{Dim,T}}) where {Dim,T} = Ring{Dim,T,typeof(vertices)}(vertices)
Ring(vertices::Tuple...) = Ring([Point(v) for v in vertices])
Ring(vertices::Point{Dim,T}...) where {Dim,T} = Ring(collect(vertices))
Ring(vertices::AbstractVector{<:Tuple}) = Ring(Point.(vertices))
Ring(vertices::AbstractVector{Point{Dim,T}}) where {Dim,T} = Ring(CircularVector(vertices))

nvertices(r::Ring) = length(r.vertices)

==(r₁::Ring, r₂::Ring) = r₁.vertices == r₂.vertices

function Base.isapprox(r₁::Ring, r₂::Ring; kwargs...)
  nvertices(r₁) ≠ nvertices(r₂) && return false
  all(isapprox(v₁, v₂; kwargs...) for (v₁, v₂) in zip(r₁.vertices, r₂.vertices))
end

Base.close(r::Ring) = r

# call `open` again to avoid issues in case of nested CircularVector
Base.open(r::Ring) = open(Rope(parent(r.vertices)))

# do not change which vertex comes first for closed chains
Base.reverse!(r::Ring) = (reverse!(@view r.vertices[(begin + 1):end]); r)

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{<:Ring{Dim,T}}) where {Dim,T}
  v = rand(rng, Point{Dim,T}, rand(3:50))
  while first(v) == last(v)
    v = rand(rng, Point{Dim,T}, rand(3:50))
  end
  Ring(v)
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

"""
    innerangles(ring)

Return inner angles of the `ring`. Inner
angles are always positive, and unlike
`angles` they can be greater than `π`.
"""
function innerangles(r::Ring{2,T}) where {T}
  # correct sign of angles in case orientation is CW
  θs = orientation(r) == CW ? -angles(r) : angles(r)
  [θ > 0 ? 2 * T(π) - θ : -θ for θ in θs]
end

innerangles(r::Ring{3}) = innerangles(Ring(proj2D(vertices(r))))
