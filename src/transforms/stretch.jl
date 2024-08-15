# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Stretch(s₁, s₂, ...)

Stretch geometry or domain outwards with
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

isrevertible(::Type{<:Stretch}) = true

isinvertible(::Type{<:Stretch}) = true

inverse(t::Stretch) = Stretch(1 ./ t.factors)

→(t₁::Stretch, t₂::Stretch) = Stretch(t₁.factors .* t₂.factors)

function apply(t::Stretch, g::GeometryOrDomain)
  p = _stretch(t, g)
  n, c = apply(p, g)
  n, (p, c)
end

revert(::Stretch, g::GeometryOrDomain, c) = revert(c[1], g, c[2])

reapply(::Stretch, g::GeometryOrDomain, c) = reapply(c[1], g, c[2])

# --------------
# SPECIAL CASES
# --------------

apply(t::Stretch, v::Vec) = apply(Scale(t.factors), v)

revert(t::Stretch, v::Vec, c) = revert(Scale(t.factors), v, c)

reapply(t::Stretch, v::Vec, c) = reapply(Scale(t.factors), v, c)

# -----------------
# HELPER FUNCTIONS
# -----------------

function _stretch(t, g)
  o = to(_origin(g))
  Translate(-o...) → Scale(t.factors) → Translate(o...)
end

_origin(g) = centroid(g)
_origin(p::Plane) = p(0, 0)
