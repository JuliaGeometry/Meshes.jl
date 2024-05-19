# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ring(p1, p2, ..., pn)

A closed polygonal chain from a sequence of points `p1`, `p2`, ..., `pn`.

See also [`Chain`](@ref) and [`Rope`](@ref).
"""
struct Ring{Dim,P<:Point{Dim},V<:CircularVector{P}} <: Chain{Dim,P}
  vertices::V

  function Ring{Dim,P,V}(vertices) where {Dim,P<:Point{Dim},V<:CircularVector{P}}
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

Ring(vertices::CircularVector{P}) where {Dim,P<:Point{Dim}} = Ring{Dim,P,typeof(vertices)}(vertices)
Ring(vertices::Tuple...) = Ring([Point(v) for v in vertices])
Ring(vertices::P...) where {P<:Point} = Ring(collect(vertices))
Ring(vertices::AbstractVector{<:Tuple}) = Ring(Point.(vertices))
Ring(vertices::AbstractVector{<:Point}) = Ring(CircularVector(vertices))

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

function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{<:Ring{Dim}}) where {Dim}
  v = [rand(rng, Point{Dim}) for _ in 1:rand(rng, 3:50)]
  while first(v) == last(v)
    v = [rand(rng, Point{Dim}) for _ in 1:rand(rng, 3:50)]
  end
  Ring(v)
end

"""
    innerangles(ring)

Return inner angles of the `ring`. Inner
angles are always positive, and unlike
`angles` they can be greater than `π`.
"""
function innerangles(r::Ring{2})
  # correct sign of angles in case orientation is CW
  θs = orientation(r) == CW ? -angles(r) : angles(r)
  [θ > 0 ? 2 * oftype(θ, π) - θ : -θ for θ in θs]
end

innerangles(r::Ring{3}) = innerangles(Ring(proj2D(vertices(r))))
