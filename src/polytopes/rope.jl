# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Rope(p1, p2, ..., pn)

An open polygonal chain from a sequence of points `p1`, `p2`, ..., `pn`.

See also [`Chain`](@ref) and [`Ring`](@ref).
"""
struct Rope{Dim,P<:Point{Dim},V<:AbstractVector{P}} <: Chain{Dim,P}
  vertices::V
end

Rope(vertices::Tuple...) = Rope([Point(v) for v in vertices])
Rope(vertices::P...) where {P<:Point} = Rope(collect(vertices))
Rope(vertices::AbstractVector{<:Tuple}) = Rope(Point.(vertices))

nvertices(r::Rope) = length(r.vertices)

==(r₁::Rope, r₂::Rope) = r₁.vertices == r₂.vertices

function Base.isapprox(r₁::Rope, r₂::Rope; kwargs...)
  nvertices(r₁) ≠ nvertices(r₂) && return false
  all(isapprox(v₁, v₂; kwargs...) for (v₁, v₂) in zip(r₁.vertices, r₂.vertices))
end

Base.close(r::Rope) = Ring(r.vertices)

Base.open(r::Rope) = r

Base.reverse!(r::Rope) = (reverse!(r.vertices); r)

Random.rand(rng::Random.AbstractRNG, ::Type{Rope{Dim}}) where {Dim} = Rope(rand(rng, Point{Dim}, rand(rng, 2:50)))
