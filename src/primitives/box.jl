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
struct Box{Dim,P<:Point{Dim}} <: Primitive{Dim}
  min::P
  max::P

  function Box{P}(min, max) where {P<:Point}
    @assert min ⪯ max "`min` must be less than or equal to `max`"
    new(min, max)
  end
end

Box(min::P, max::P) where {P<:Point} = Box{P}(min, max)

Box(min::Tuple, max::Tuple) = Box(Point(min), Point(max))

paramdim(::Type{<:Box{Dim}}) where {Dim} = Dim

Base.minimum(b::Box) = b.min

Base.maximum(b::Box) = b.max

Base.extrema(b::Box) = b.min, b.max

center(b::Box) = Point((coordinates(b.max) + coordinates(b.min)) / 2)

diagonal(b::Box) = norm(b.max - b.min)

sides(b::Box) = Tuple(b.max - b.min)

Base.isapprox(b₁::Box, b₂::Box) = b₁.min ≈ b₂.min && b₁.max ≈ b₂.max

# TODO
# function (b::Box{Dim,T})(uv...) where {Dim,T}
#   if !all(x -> zero(T) ≤ x ≤ one(T), uv)
#     throw(DomainError(uv, "b(u, v, ...) is not defined for u, v, ... outside [0, 1]ⁿ."))
#   end
#   b.min + uv .* (b.max - b.min)
# end

# TODO
# function Random.rand(rng::Random.AbstractRNG, ::Random.SamplerType{Box{Dim,T}}) where {Dim,T}
#   min = rand(rng, Point{Dim,T})
#   max = min + rand(rng, Vec{Dim,T})
#   Box(min, max)
# end
