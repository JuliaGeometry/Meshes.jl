# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Ring(p1, p2, ..., pn)

A closed polygonal chain from a sequence of points `p1`, `p2`, ..., `pn`.

See also [`Chain`](@ref) and [`Rope`](@ref).
"""
struct Ring{Dim,C<:CRS,V<:CircularVector{Point{Dim,C}}} <: Chain{Dim,C}
  vertices::V

  function Ring{Dim,C,V}(vertices) where {Dim,C<:CRS,V<:CircularVector{Point{Dim,C}}}
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

Ring(vertices::CircularVector{Point{Dim,C}}) where {Dim,C<:CRS} = Ring{Dim,C,typeof(vertices)}(vertices)
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

# call `open` again to avoid issues in case of nested CircularVector
Base.open(r::Ring) = open(Rope(parent(r.vertices)))

# do not change which vertex comes first for closed chains
Base.reverse!(r::Ring) = (reverse!(@view r.vertices[(begin + 1):end]); r)

function Random.rand(rng::Random.AbstractRNG, ::Type{Ring{Dim}}) where {Dim}
  v = rand(rng, Point{Dim}, rand(rng, 3:50))
  while first(v) == last(v)
    v = rand(rng, Point{Dim}, rand(rng, 3:50))
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
