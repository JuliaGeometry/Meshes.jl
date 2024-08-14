# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Box(min, max)

A (geodesic) box with `min` and `max` points on a given manifold.

## Examples

Construct a 3D box using points with Cartesian coordinates:

```julia
Box((0, 0, 0), (1, 1, 1))
```

Likewise, construct a 2D box on the plane:

```julia
Box((0, 0), (1, 1))
```

Construct a geodesic box on the ellipsoid:

```julia
Box(Point(LatLon(0, 0)), Point(LatLon(1, 1)))
```
"""
struct Box{M<:Manifold,C<:CRS} <: Primitive{M,C}
  min::Point{M,C}
  max::Point{M,C}

  function Box{M,C}(min, max) where {M<:Manifold,C<:CRS}
    assertion(min ⪯ max, "can only construct box with min ⪯ max")
    new(min, max)
  end
end

Box(min::Point{M,C}, max::Point{M,C}) where {M<:Manifold,C<:CRS} = Box{M,C}(min, max)

Box(min::Tuple, max::Tuple) = Box(Point(min), Point(max))

paramdim(::Type{<:Box{𝔼{Dim}}}) where {Dim} = Dim

paramdim(::Type{<:Box{🌐}}) = 2

Base.minimum(b::Box) = b.min

Base.maximum(b::Box) = b.max

Base.extrema(b::Box) = b.min, b.max

diagonal(b::Box) = norm(b.max - b.min)

sides(b::Box) = Tuple(b.max - b.min)

==(b₁::Box, b₂::Box) = b₁.min == b₂.min && b₁.max == b₂.max

Base.isapprox(b₁::Box, b₂::Box; atol=atol(lentype(b₁)), kwargs...) =
  isapprox(b₁.min, b₂.min; atol, kwargs...) && isapprox(b₁.max, b₂.max; atol, kwargs...)

function (b::Box{𝔼})(uv...)
  if !all(x -> 0 ≤ x ≤ 1, uv)
    throw(DomainError(uv, "b(u, v, ...) is not defined for u, v, ... outside [0, 1]ⁿ."))
  end
  b.min + uv .* (b.max - b.min)
end
