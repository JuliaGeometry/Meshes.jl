# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Box(min, max)

An axis-aligned box with `min` and `max` corners.
See <https://en.wikipedia.org/wiki/Hyperrectangle>.

## Examples

```julia
Box(Point(0, 0, 0), Point(1, 1, 1))
Box((0, 0), (1, 1))
```
"""
struct Box{Dim,C<:CRS} <: Primitive{Dim,C}
  min::Point{Dim,C}
  max::Point{Dim,C}

  function Box{Dim,C}(min, max) where {Dim,C<:CRS}
    assertion(min ⪯ max, "`min` must be less than or equal to `max`")
    new(min, max)
  end
end

Box(min::Point{Dim,C}, max::Point{Dim,C}) where {Dim,C<:CRS} = Box{Dim,C}(min, max)

Box(min::Tuple, max::Tuple) = Box(Point(min), Point(max))

paramdim(::Type{<:Box{Dim}}) where {Dim} = Dim

Base.minimum(b::Box) = b.min

Base.maximum(b::Box) = b.max

Base.extrema(b::Box) = b.min, b.max

center(b::Box) = withdatum(b, (to(b.max) + to(b.min)) / 2)

diagonal(b::Box) = norm(b.max - b.min)

sides(b::Box) = Tuple(b.max - b.min)

==(b₁::Box, b₂::Box) = b₁.min == b₂.min && b₁.max == b₂.max

Base.isapprox(b₁::Box, b₂::Box; atol=atol(lentype(b₁)), kwargs...) =
  isapprox(b₁.min, b₂.min; atol, kwargs...) && isapprox(b₁.max, b₂.max; atol, kwargs...)

function (b::Box)(uv...)
  if !all(x -> 0 ≤ x ≤ 1, uv)
    throw(DomainError(uv, "b(u, v, ...) is not defined for u, v, ... outside [0, 1]ⁿ."))
  end
  b.min + uv .* (b.max - b.min)
end

function Random.rand(rng::Random.AbstractRNG, ::Type{Box{Dim}}) where {Dim}
  min = rand(rng, Point{Dim})
  max = min + rand(rng, Vec{Dim,Met{Float64}})
  Box(min, max)
end
