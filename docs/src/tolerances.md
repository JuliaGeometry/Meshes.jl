# Tolerances

The absolute tolerance used for floating point arithmetic is hard-coded in
the project to `1e-10` for `Float64` and to `1f-5` for `Float32`. You can use
[ScopedValues.jl](https://github.com/vchuravy/ScopedValues.jl) to customize
these tolerance values in specific computations:

```julia
using Meshes
using ScopedValues

with(Meshes.ATOL64 => 1e-9, Meshes.ATOL32 => 1f-4) do
  # do your computations with custom tolerances
end
```
