# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Shadow(plane)

TODO

## Examples

```julia
Shadow(:xy)
Shadow("xz")
```
"""
struct Shadow{Dim} <: CoordinateTransform
  plane::NTuple{Dim,Symbol}
end

Shadow(plane::AbstractString) = Shadow(Tuple(Symbol.(sort(split(plane, "")))))
Shadow(plane::Symbol) = Shadow(string(plane))

parameters(t::Shadow) = (; plane=t.plane)

_index(::Val{:x}) = 1
_index(::Val{:y}) = 2
_index(::Val{:z}) = 3

applycoord(t::Shadow{Dim}, v::Vec) where {Dim} = Vec(ntuple(i -> v[_index(Val(t.plane[i]))], Dim))

function applycoord(t::Shadow{Dim}, p::Point) where {Dim}
  c = convert(Cartesian, coords(p))
  Point(Cartesian{datum(c)}(ntuple(i -> getproperty(c, t.plane[i]), Dim)))
end

# --------------
# SPECIAL CASES
# --------------

applycoord(t::Shadow, g::RectilinearGrid) = applycoord(t, convert(SimpleMesh, g))

applycoord(t::Shadow, g::StructuredGrid) = applycoord(t, convert(SimpleMesh, g))
