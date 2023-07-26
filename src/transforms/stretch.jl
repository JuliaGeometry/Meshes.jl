# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Stretch(s₁, s₂, ...)

Scale coordinates of geometry or mesh by
given strictly positive factors `s₁, s₂, ...`.

## Examples

```julia
Stretch(1.0, 2.0, 3.0)
```
"""
struct Stretch{Dim,T} <: CoordinateTransform
  factors::NTuple{Dim,T}

  function Stretch{Dim,T}(factors) where {Dim,T}
    any(≤(0), factors) && throw(ArgumentError("Scaling factors must be positive."))
    new(factors)
  end
end

Stretch(factors::NTuple{Dim,T}) where {Dim,T} = Stretch{Dim,T}(factors)

Stretch(factors...) = Stretch(factors)

isrevertible(::Type{<:Stretch}) = true

isinvertible(::Type{<:Stretch}) = true

Base.inv(t::Stretch) = Stretch(1 ./ t.factors)

applycoord(t::Stretch, v::Vec) = t.factors .* v

# --------------
# SPECIAL CASES
# --------------

function applycoord(t::Stretch, g::CartesianGrid)
  dims = size(g)
  orig = applycoord(t, minimum(g))
  spac = t.factors .* spacing(g)
  offs = offset(g)
  CartesianGrid(dims, orig, spac, offs)
end