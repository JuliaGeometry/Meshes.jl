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
struct Stretch{Dim,T} <: PointwiseGeometricTransform
  factors::NTuple{Dim,T}

  function Stretch{Dim,T}(factors) where {Dim,T}
    if any(≤(0), factors)
      throw(ArgumentError("Scaling factors must be positive."))
    end
    new(factors)
  end
end

Stretch(factors::NTuple{Dim,T}) where {Dim,T} = Stretch{Dim,T}(factors)

Stretch(factors...) = Stretch(factors)

Base.inv(t::Stretch) = Stretch(1 ./ t.factors)

isrevertible(::Type{<:Stretch}) = true

_apply(t::Stretch, v::Vec) = t.factors .* v
