# Tolerances

## Absolute tolerance for floating point comparisons

The absolute tolerance used for floating point comparisons is hard-coded in
the project to `eps(T) ^ (3 // 4)` where `T` is either `Float64` or `Float32`.
This formula preserves roughly 3/4 of the significant digits of the type `T`.
You can use [ScopedValues.jl](https://github.com/vchuravy/ScopedValues.jl) to
customize these tolerance values in specific computations:

```julia
using Meshes
using ScopedValues

with(Meshes.ATOL64 => 1e-9, Meshes.ATOL32 => 1f-4) do
  # do your computations with custom tolerances
end
```

## Maximum length for default refinement of meshes on `🌐`

The default discretization of geometries on the `🌐` manifold includes an additional
refinement step to make sure that all segments of the mesh are shorter than or equal
to a maximum length of `500km`. This value can be customized in specific computations
or visualizations:

```julia
using Meshes
using ScopedValues

with(Meshes.MAXLEN => 100km) do
  # do your computations with custom maximum length
end
```
