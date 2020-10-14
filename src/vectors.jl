"""
    Vec{N,T}

A vector in `N`-dimensional space with coordinates of type
`T` representing a direction with magnitude. A vector can
be obtained by subtracting two [`Point`](@ref) objects:

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
    vunit(VecType, i)

Return the `i`-th vector in the Euclidean basis
as a vector with type `VecType`.

## Example

```julia
vunit(Vec3, 1) == [1.0, 0.0, 0.0]
vunit(Vec3, 2) == [0.0, 1.0, 0.0]
vunit(Vec3, 3) == [0.0, 0.0, 1.0]
```
"""
function vunit(VecType::Type, i::Integer)
    N = length(VecType)
    T = eltype(VecType)
    VecType(ntuple(j->ifelse(i==j, one(T), zero(T)), N))
end

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
