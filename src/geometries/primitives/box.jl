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

Base.extrema(b::Box) = minimum(b), maximum(b)

diagonal(b::Box{<:𝔼}) = norm(maximum(b) - minimum(b))

sides(b::Box{<:𝔼}) = Tuple(maximum(b) - minimum(b))

==(b₁::Box, b₂::Box) = minimum(b₁) == minimum(b₂) && maximum(b₁) == maximum(b₂)

Base.isapprox(b₁::Box, b₂::Box; atol=atol(lentype(b₁)), kwargs...) =
  isapprox(minimum(b₁), minimum(b₂); atol, kwargs...) && isapprox(maximum(b₁), maximum(b₂); atol, kwargs...)

(b::Box{<:𝔼})(uvw...) = minimum(b) + promote(uvw...) .* (maximum(b) - minimum(b))
