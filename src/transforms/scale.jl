# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Scale(s₁, s₂, ...)

Scale geometry or domain with
strictly positive scaling factors
`s₁, s₂, ...`.

## Examples

```julia
Scale(1.0, 2.0, 3.0)
```
"""
struct Scale{Dim,T} <: CoordinateTransform
  factors::NTuple{Dim,T}

  function Scale{Dim,T}(factors) where {Dim,T}
    any(≤(0), factors) && throw(ArgumentError("Scaling factors must be positive."))
    new(factors)
  end
end

Scale(factors::NTuple{Dim,T}) where {Dim,T} = Scale{Dim,T}(factors)

Scale(factors...) = Scale(factors)

parameters(t::Scale) = (; factors=t.factors)

isaffine(::Type{<:Scale}) = true

isrevertible(::Type{<:Scale}) = true

isinvertible(::Type{<:Scale}) = true

inverse(t::Scale) = Scale(1 ./ t.factors)

applycoord(t::Scale, v::Vec) = t.factors .* v

# --------------
# SPECIAL CASES
# --------------

applycoord(t::Scale, b::Ball) = applycoord(t, discretize(b))

applycoord(t::Scale{Dim}, s::Sphere) where {Dim} = _applycoord(t, s, Val(Dim), Val(embeddim(s)))

_applycoord(t::Scale, s::Sphere, ::Val, ::Val) = applycoord(t, discretize(s))

_applycoord(t::Scale, s::Sphere, ::Val{1}, ::Val) = Sphere(applycoord(t, center(s)), t.factors[1] * radius(s))

_applycoord(t::Scale, s::Sphere, ::Val{3}, ::Val{3}) = Ellipsoid(t.factors .* radius(s), applycoord(t, center(s)))

applycoord(t::Scale, e::Ellipsoid) = applycoord(t, discretize(e))

applycoord(t::Scale, d::Disk) = applycoord(t, discretize(d))

applycoord(t::Scale, c::Circle) = applycoord(t, discretize(c))

applycoord(t::Scale, c::Cylinder) = applycoord(t, discretize(c))

applycoord(t::Scale, c::CylinderSurface) = applycoord(t, discretize(c))

applycoord(t::Scale, p::ParaboloidSurface) = applycoord(t, discretize(p))

applycoord(t::Scale, tr::Torus) = applycoord(t, discretize(tr))

function applycoord(t::Scale, g::CartesianGrid)
  dims = size(g)
  orig = applycoord(t, minimum(g))
  spac = t.factors .* spacing(g)
  offs = offset(g)
  CartesianGrid(dims, orig, spac, offs)
end

applycoord(t::Scale{Dim}, g::RectilinearGrid{Datum,Dim}) where {Datum,Dim} =
  RectilinearGrid{Datum}(ntuple(i -> t.factors[i] * xyz(g)[i], Dim))

applycoord(t::Scale{Dim}, g::StructuredGrid{Datum,Dim}) where {Datum,Dim} =
  StructuredGrid{Datum}(ntuple(i -> t.factors[i] * XYZ(g)[i], Dim))
