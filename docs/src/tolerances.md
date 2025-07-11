# Tolerances

The absolute tolerance used for floating point comparisons is hard-coded in
the project to `eps(T) ^ (3 // 4)` where `T` is either `Float64` or `Float32`.
You can use [ScopedValues.jl](https://github.com/vchuravy/ScopedValues.jl) to
customize these tolerance values in specific computations:

```julia
using Meshes
using ScopedValues

with(Meshes.ATOL64 => 1e-9, Meshes.ATOL32 => 1f-4) do
  # do your computations with custom tolerances
end
```
