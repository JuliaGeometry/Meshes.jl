# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ring(p1, p2, ..., pn)

A closed polygonal chain from a sequence of points `p1`, `p2`, ..., `pn`.

See also [`Chain`](@ref) and [`Rope`](@ref).
"""
struct Ring{M<:Manifold,C<:CRS,V<:CircularVector{Point{M,C}}} <: Chain{M,C}
  vertices::V

  function Ring{M,C,V}(vertices) where {M<:Manifold,C<:CRS,V<:CircularVector{Point{M,C}}}
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

Ring(vertices::CircularVector{Point{M,C}}) where {M<:Manifold,C<:CRS} = Ring{M,C,typeof(vertices)}(vertices)
Ring(vertices::Tuple...) = Ring([Point(v) for v in vertices])
Ring(vertices::P...) where {P<:Point} = Ring(collect(vertices))
Ring(vertices::AbstractVector{<:Tuple}) = Ring(Point.(vertices))
Ring(vertices::AbstractVector{<:Point}) = Ring(CircularVector(vertices))

nvertices(r::Ring) = length(r.vertices)

==(r₁::Ring, r₂::Ring) = r₁.vertices == r₂.vertices

"""
    ≗(ring₁, ring₂)

Tells whether or not the `ring₁` and `ring₂`
are equal up to circular shifts.
"""
function ≗(r₁::Ring, r₂::Ring)
  n = length(r₁.vertices)
  i = findfirst(==(first(r₁.vertices)), r₂.vertices)
  isnothing(i) && return false
  r₁.vertices == r₂.vertices[i:(i + n - 1)]
end

Base.isapprox(r₁::Ring, r₂::Ring; atol=atol(lentype(r₁)), kwargs...) =
  nvertices(r₁) == nvertices(r₂) && all(isapprox(v₁, v₂; atol, kwargs...) for (v₁, v₂) in zip(r₁.vertices, r₂.vertices))

Base.close(r::Ring) = r

# call `open` again to avoid issues with nested CircularVector
Base.open(r::Ring) = open(Rope(parent(r.vertices)))

# do not change which vertex comes first in reverse order
Base.reverse!(r::Ring) = (reverse!(@view r.vertices[(begin + 1):end]); r)

"""
    innerangles(ring)

Return inner angles of the `ring`. Inner
angles are always positive, and unlike
`angles` they can be greater than `π`.
"""
function innerangles(r::Ring)
  # correct sign of angles in case orientation is CW
  θs = orientation(r) == CW ? -angles(r) : angles(r)
  [θ > 0 ? 2 * oftype(θ, π) - θ : -θ for θ in θs]
end
