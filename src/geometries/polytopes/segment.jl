# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Segment(p1, p2)

An oriented, geodesic line segment from point `p1` to point `p2`.

## Examples

```julia
# straight segment in Euclidean space
Segment((0, 0), (1, 1))

# geodesic segment in Earth's surface
Segment(Point(LatLon(0, 0)), Point(LatLon(45, 90)))
```

See also [`Rope`](@ref), [`Ring`](@ref), [`Line`](@ref).
"""
@polytope Segment 1 2

nvertices(::Type{<:Segment}) = 2

Base.minimum(s::Segment) = s.vertices[1]

Base.maximum(s::Segment) = s.vertices[2]

Base.extrema(s::Segment) = s.vertices[1], s.vertices[2]

center(s::Segment{<:𝔼}) = s(1 // 2)

==(s₁::Segment, s₂::Segment) = s₁.vertices == s₂.vertices

Base.isapprox(s₁::Segment, s₂::Segment; atol=atol(lentype(s₁)), kwargs...) =
  all(isapprox(v₁, v₂; atol, kwargs...) for (v₁, v₂) in zip(s₁.vertices, s₂.vertices))

function (s::Segment{<:𝔼})(t)
  a, b = s.vertices
  a + t * (b - a)
end

function (s::Segment{🌐})(t)
  a, b = s.vertices
  d = GeodesicDistance()
  ϕ = geodesicbwd(a, b)
  geodesicfwd(a, ϕ, t * d(a, b))
end

Base.reverse(s::Segment) = Segment(reverse(s.vertices))
