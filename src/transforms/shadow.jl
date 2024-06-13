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
  plane::NTuple{Dim,Symbol}
end

Shadow(plane::AbstractString) = Shadow(Tuple(Symbol.(sort(split(plane, "")))))
Shadow(plane::Symbol) = Shadow(string(plane))

parameters(t::Shadow) = (; plane=t.plane)

apply(t::Shadow, v::Vec) = _shadow(t, v), nothing

apply(t::Shadow, g::GeometryOrDomain) = _shadow(t, g), nothing

_index(::Val{:x}) = 1
_index(::Val{:y}) = 2
_index(::Val{:z}) = 3

function _shadow(t::Shadow{Dim}, g) where {Dim}
  inds = SVector(ntuple(i -> _index(Val(t.plane[i])), Dim))
  _shadowinds(g, inds)
end

_shadowinds(v::Vec, inds) = v[inds]

_shadowinds(p::Point, inds) = withdatum(p, _shadowinds(to(p), inds))

function _shadowinds(g::CartesianGrid, inds)
  sz = size(g)[inds]
  or = _shadowinds(minimum(g), inds)
  sp = spacing(g)[inds]
  of = offset(g)[inds]
  CartesianGrid(sz, or, sp, of)
end

_shadowinds(g::RectilinearGrid, inds) = RectilinearGrid{datum(crs(g))}(xyz(g)[inds])

_dropdims(X, inds) = X[ntuple(i -> ifelse(i ∈ inds, :, 1), ndims(X))...]

_shadowinds(g::StructuredGrid, inds) = StructuredGrid{datum(crs(g))}(map(X -> _dropdims(X, inds), XYZ(g)[inds]))

# apply shadow transform recursively
@generated function _shadowinds(g::G, inds) where {G<:GeometryOrDomain}
  ctor = constructor(G)
  names = fieldnames(G)
  exprs = (:(_shadowinds(g.$name, inds)) for name in names)
  :($ctor($(exprs...)))
end

# stop recursion at non-geometric types
_shadowinds(x, _) = x

# special treatment for lists of geometries
_shadowinds(g::NTuple{<:Any,<:Geometry}, inds) = map(gᵢ -> _shadowinds(gᵢ, inds), g)
_shadowinds(g::AbstractVector{<:Geometry}, inds) = tcollect(_shadowinds(gᵢ, inds) for gᵢ in g)
_shadowinds(g::CircularVector{<:Geometry}, inds) = CircularVector(tcollect(_shadowinds(gᵢ, inds) for gᵢ in g))
