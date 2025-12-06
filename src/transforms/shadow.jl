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

apply(t::Shadow, v::Vec) = v[_sort(t.dims)], nothing

function apply(t::Shadow, g::GeometryOrDomain)
  dims = _sort(t.dims)
  m = Morphological() do coords
    cart = convert(Cartesian, coords)
    vals = CoordRefSystems.values(cart)
    Cartesian{datum(coords)}(vals[dims])
  end
  apply(m, g)
end

# --------------
# SPECIAL CASES
# --------------

apply(t::Shadow, b::Box{<:𝔼}) = Box(t(minimum(b)), t(maximum(b))), nothing

apply(::Shadow, ::Plane) = throw(ArgumentError("Shadow transform doesn't yet support planes"))

function apply(t::Shadow, g::CartesianGrid)
  dims = _sort(t.dims)
  or = t(minimum(g))
  sp = spacing(g)[dims]
  sz = size(g)[dims]
  ip = isperiodic(g)[dims]
  tp = GridTopology(sz, ip)
  CartesianGrid(or, sp, tp), nothing
end

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
