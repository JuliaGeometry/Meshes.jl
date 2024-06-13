# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Shadow(plane)

Project the geometry or domain onto the given `plane`,
producing a "shadow" of the original object.

## Examples

```julia
Shadow(:xy)
Shadow("xz")
```
"""
struct Shadow{Dim} <: GeometricTransform
  inds::SVector{Dim,Int}
  Shadow(inds::SVector{Dim,Int}) where {Dim} = new{Dim}(sort(inds))
end

Shadow(inds::NTuple{Dim,Int}) where {Dim} = Shadow(SVector(inds))

Shadow(inds::Int...) = Shadow(inds)

Shadow(plane::AbstractString) = Shadow(Tuple(_index(Val(Symbol(s))) for s in eachsplit(plane, "")))

Shadow(plane::Symbol) = Shadow(string(plane))

parameters(t::Shadow) = (; inds=t.inds)

apply(t::Shadow, v::Vec) = _shadow(v, t.inds), nothing

apply(t::Shadow, g::GeometryOrDomain) = _shadow(g, t.inds), nothing

_index(::Val{:x}) = 1
_index(::Val{:y}) = 2
_index(::Val{:z}) = 3

_shadow(v::Vec, inds) = v[inds]

_shadow(p::Point, inds) = withdatum(p, to(p)[inds])

function _shadow(g::CartesianGrid, inds)
  sz = size(g)[inds]
  or = _shadow(minimum(g), inds)
  sp = spacing(g)[inds]
  of = offset(g)[inds]
  CartesianGrid(sz, or, sp, of)
end

_shadow(g::RectilinearGrid, inds) = RectilinearGrid{datum(crs(g))}(xyz(g)[inds])

function _shadow(g::StructuredGrid, inds)
  ndims = length(size(g))
  slices = ntuple(i -> ifelse(i ∈ inds, :, 1), ndims)
  StructuredGrid{datum(crs(g))}(map(X -> X[slices...], XYZ(g)[inds]))
end

# apply shadow transform recursively
@generated function _shadow(g::G, inds) where {G<:GeometryOrDomain}
  ctor = constructor(G)
  names = fieldnames(G)
  exprs = (:(_shadow(g.$name, inds)) for name in names)
  :($ctor($(exprs...)))
end

# stop recursion at non-geometric types
_shadow(x, _) = x

# special treatment for lists of geometries
_shadow(g::NTuple{<:Any,<:Geometry}, inds) = map(gᵢ -> _shadow(gᵢ, inds), g)
_shadow(g::AbstractVector{<:Geometry}, inds) = tcollect(_shadow(gᵢ, inds) for gᵢ in g)
_shadow(g::CircularVector{<:Geometry}, inds) = CircularVector(tcollect(_shadow(gᵢ, inds) for gᵢ in g))
