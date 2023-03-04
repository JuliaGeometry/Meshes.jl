# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Stretch(s₁, s₂, ...)

Scale coordinates of geometry or mesh by
given factors `s₁, s₂, ...`.

## Examples

```julia
Stretch(1.0, 2.0, 3.0)
```
"""
struct Stretch{Dim,T} <: StatelessGeometricTransform
  factors::NTuple{Dim,T}

  function Stretch{Dim,T}(factors) where {Dim,T}
    if any(≤(0), factors)
      throw(ArgumentError("Scaling factors must be positive."))
    end
    return new(factors)
  end
end

Stretch(factors::NTuple{Dim,T}) where {Dim,T} = Stretch{Dim,T}(factors)

Stretch(factors...) = Stretch(factors)

isrevertible(::Type{<:Stretch}) = true

preprocess(transform::Stretch, object) = transform.factors

function applypoint(::Stretch, points, prep)
  s = prep
  newpoints = [Point(s .* coordinates(p)) for p in points]
  return newpoints, prep
end

function revertpoint(::Stretch, newpoints, cache)
  s = cache
  return [Point((1 ./ s) .* coordinates(p)) for p in newpoints]
end
