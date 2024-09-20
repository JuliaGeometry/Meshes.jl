# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

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

# --------------
# SPECIAL CASES
# --------------

apply(::Shadow, ::Plane) = throw(ArgumentError("Shadow transform doesn't yet support planes"))

apply(t::Shadow, e::Ellipsoid) = TransformedGeometry(e, t), nothing

apply(t::Shadow, d::Disk) = TransformedGeometry(d, t), nothing

apply(t::Shadow, c::Circle) = TransformedGeometry(c, t), nothing

apply(t::Shadow, c::Cylinder) = TransformedGeometry(c, t), nothing

apply(t::Shadow, c::CylinderSurface) = TransformedGeometry(c, t), nothing

apply(t::Shadow, c::Cone) = TransformedGeometry(c, t), nothing

apply(t::Shadow, c::ConeSurface) = TransformedGeometry(c, t), nothing

apply(t::Shadow, f::Frustum) = TransformedGeometry(f, t), nothing

apply(t::Shadow, f::FrustumSurface) = TransformedGeometry(f, t), nothing

apply(t::Shadow, p::ParaboloidSurface) = TransformedGeometry(p, t), nothing

apply(t::Shadow, tr::Torus) = TransformedGeometry(tr, t), nothing

apply(t::Shadow, ct::CylindricalTrajectory) = apply(t, GeometrySet(collect(ct))), nothing

apply(t::Shadow, g::CartesianGrid) = _shadow(g, _sort(t.dims)), nothing

apply(t::Shadow, g::RegularGrid) = TransformedGrid(g, t), nothing

apply(t::Shadow, g::RectilinearGrid) = TransformedGrid(g, t), nothing

apply(t::Shadow, g::StructuredGrid) = TransformedGrid(g, t), nothing

# -----------------
# HELPER FUNCTIONS
# -----------------

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

_sort(dims) = sort(SVector(dims))

_shadow(v::Vec, dims) = v[dims]

function _shadow(p::Point, dims)
  v = _shadow(to(p), dims)
  c = Cartesian{datum(crs(p))}(v...)
  Point(c)
end

function _shadow(g::CartesianGrid, dims)
  sz = size(g)[dims]
  or = _shadow(minimum(g), dims)
  sp = spacing(g)[dims]
  of = offset(g)[dims]
  CartesianGrid(sz, or, sp, of)
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
_shadow(g::AbstractVector{<:Geometry}, dims) = [_shadow(gᵢ, dims) for gᵢ in g]
_shadow(g::CircularVector{<:Geometry}, dims) = CircularVector([_shadow(gᵢ, dims) for gᵢ in g])
