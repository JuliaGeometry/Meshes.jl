# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function _index(d)
  if d == 'x'
    1
  elseif d == 'y'
    2
  elseif d == 'z'
    3
  else
    throw(ArgumentError("'$d' isn't a valid dimension name"))
  end
end

"""
    Shadow(dims)

Project the geometry or domain onto the given `dims`,
producing a "shadow" of the original object.

## Examples

```julia
Shadow(:xy)
Shadow("xz")
Shadow(1, 2)
Shadow((1, 3))
```
"""
struct Shadow{Dim} <: GeometricTransform
  dims::Dims{Dim}
end

Shadow(dims::Int...) = Shadow(dims)

Shadow(dims::AbstractString) = Shadow(Dims(_index(d) for d in dims))

Shadow(dims::Symbol) = Shadow(string(dims))

parameters(t::Shadow) = (; dims=t.dims)

apply(t::Shadow, v::Vec) = _shadow(v, _sort(t.dims)), nothing

apply(t::Shadow, g::GeometryOrDomain) = _shadow(g, _sort(t.dims)), nothing

_sort(dims) = sort(SVector(dims))

_shadow(v::Vec, dims) = v[dims]

_shadow(p::Point, dims) = withdatum(p, to(p)[dims])

function _shadow(g::CartesianGrid, dims)
  sz = size(g)[dims]
  or = _shadow(minimum(g), dims)
  sp = spacing(g)[dims]
  of = offset(g)[dims]
  CartesianGrid(sz, or, sp, of)
end

_shadow(g::RectilinearGrid, dims) = RectilinearGrid{datum(crs(g))}(xyz(g)[dims])

function _shadow(g::StructuredGrid, dims)
  ndims = length(size(g))
  inds = ntuple(i -> ifelse(i ∈ dims, :, 1), ndims)
  StructuredGrid{datum(crs(g))}(map(X -> X[inds...], XYZ(g)[dims]))
end

# apply shadow transform recursively
@generated function _shadow(g::G, dims) where {G<:GeometryOrDomain}
  ctor = constructor(G)
  names = fieldnames(G)
  exprs = (:(_shadow(g.$name, dims)) for name in names)
  :($ctor($(exprs...)))
end

# stop recursion at non-geometric types
_shadow(x, _) = x

# special treatment for lists of geometries
_shadow(g::NTuple{<:Any,<:Geometry}, dims) = map(gᵢ -> _shadow(gᵢ, dims), g)
_shadow(g::AbstractVector{<:Geometry}, dims) = tcollect(_shadow(gᵢ, dims) for gᵢ in g)
_shadow(g::CircularVector{<:Geometry}, dims) = CircularVector(tcollect(_shadow(gᵢ, dims) for gᵢ in g))
