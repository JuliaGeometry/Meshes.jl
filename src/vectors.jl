# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    Vec{Dim,T}

A vector in `Dim`-dimensional space with coordinates of type `T`
representing a direction with magnitude. A vector can be obtained
by subtracting two [`Point`](@ref) objects:

## Example

```julia
A = Point(0.0, 0.0)
B = Point(1.0, 0.0)
v = B - A
```

### Notes

- A `Vec` is a `SVector` from StaticArrays.jl
- Type aliases are `Vec2`, `Vec3`, `Vec2f`, `Vec3f`
"""
const Vec = SVector

# type aliases for convenience
const Vec2  = Vec{2,Float64}
const Vec3  = Vec{3,Float64}
const Vec2f = Vec{2,Float32}
const Vec3f = Vec{3,Float32}

"""
    vfill(VecType, v)

Return the vector of type `VecType` with all
coordinates equal to `v`.

## Example

```julia
vfill(Vec2, 3) == [3.0, 3.0]
vfill(Vec2f, 0) == [0.0f0, 0.0f0]
```
"""
function vfill(VecType::Type, v::Number)
  VecType(ntuple(j->v, length(VecType)))
end
