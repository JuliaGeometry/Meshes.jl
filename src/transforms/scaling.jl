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
struct ScaleCoords{T} <: GeometricTransform
  factors::Vector{T}
  function ScaleCoords{T}(factors) where {T}
    if any(≤(0), factors)
      throw(ArgumentError("Scaling factors must be positive."))
    end
    new(factors)
  end
end

ScaleCoords(factors...) = ScaleCoords{eltype(factors)}(collect(factors))

isrevertible(::Type{<:ScaleCoords}) = true

preprocess(transform::ScaleCoords, object) = transform.factors

function applypoint(::ScaleCoords, points, prep)
  s = prep
  newpoints = [Point(s .* coordinates(p)) for p in points]
  newpoints, prep
end

function revertpoint(::ScaleCoords, newpoints, cache)
  s = cache
  [Point((1 ./ s) .* coordinates(p)) for p in newpoints]
end

function reapplypoint(::ScaleCoords, points, cache)
  s = cache
  [Point(s .* coordinates(p)) for p in points]
end
