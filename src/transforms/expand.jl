# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Expand(s₁, s₂, ...)

Expand geometry or domain outwards with
strictly positive scaling factors
`s₁, s₂, ...`.

## Examples

```julia
Expand(1.0, 2.0, 3.0)
```
"""
struct Expand{Dim,T} <: CoordinateTransform
  factors::NTuple{Dim,T}

  function Expand{Dim,T}(factors) where {Dim,T}
    any(≤(0), factors) && throw(ArgumentError("Scaling factors must be positive."))
    new(factors)
  end
end

Expand(factors::NTuple{Dim,T}) where {Dim,T} = Expand{Dim,T}(factors)

Expand(factors...) = Expand(factors)

parameters(t::Expand) = (; factors=t.factors)

isrevertible(::Type{<:Expand}) = true

isinvertible(::Type{<:Expand}) = true

inverse(t::Expand) = Expand(1 ./ t.factors)

function apply(t::Expand, g::GeometryOrDomain)
  p = _expand(t, g)
  n, c = apply(p, g)
  n, (p, c)
end

revert(t::Expand, g::GeometryOrDomain, c) = revert(c[1], g, c[2])

reapply(t::Expand, g::GeometryOrDomain, c) = reapply(c[1], g, c[2])

function _expand(t, g)
  o = coordinates(_origin(g))
  Translate(-o...) → Stretch(t.factors) → Translate(o...)
end

_origin(g) = centroid(g)
_origin(p::Plane) = p(0, 0)

# --------------
# SPECIAL CASES
# --------------

apply(t::Expand, v::Vec) = apply(Stretch(t.factors), v)

revert(t::Expand, v::Vec, c) = revert(Stretch(t.factors), v, c)

reapply(t::Expand, v::Vec, c) = reapply(Stretch(t.factors), v, c)
