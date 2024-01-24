# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Stretch(s₁, s₂, ...)

Stretch geometry or domain with
strictly positive scaling factors
`s₁, s₂, ...`.

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

parameters(t::Stretch) = (; factors=t.factors)

isaffine(::Type{<:Stretch}) = true

isrevertible(::Type{<:Stretch}) = true

isinvertible(::Type{<:Stretch}) = true

inverse(t::Stretch) = Stretch(1 ./ t.factors)

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
