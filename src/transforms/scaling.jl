# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    ScaleCoords(sx, sy, sz)

Transform geometry or mesh by scaling each coordinate of the vertices by 
the corresponding factor. 

## Arguments

- `sx`: scaling factor for the first coordinate
- `sy`: scaling factor for the second coordinate
- `sz`: scaling factor for the third coordinate

## Examples

```julia
ScaleCoords(1, 2, 3)
```
"""
struct ScaleCoords{T<:Real} <: GeometricTransform
  factors::AbstractVector{T}
  function ScaleCoords{T}(factors::AbstractVector{T}) where {T<:Real}
    if any(factors .≤ 0)
      error("The scaling factors must be nonnegative.")
    end
    new(factors)
  end
end

ScaleCoords(factors...) = ScaleCoords{eltype(factors)}([factors...])

isrevertible(::Type{<:ScaleCoords}) = true

preprocess(transform::ScaleCoords, object) = transform.factors

function applypoint(::ScaleCoords, points, prep)
  scale_vector = prep
  newpoints = [Point(scale_vector .* coordinates(p)) for p in points]
  newpoints, prep
end

function revertpoint(::ScaleCoords, newpoints, cache)
  scale_vector = cache
  iscale_vector = 1 ./ scale_vector
  [Point(iscale_vector .* coordinates(p)) for p in newpoints]
end

function reapplypoint(::ScaleCoords, points, cache)
  scale_vector = cache
  [Point(scale_vector .* coordinates(p)) for p in points]
end
