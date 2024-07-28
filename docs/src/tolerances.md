# Changing default absolute tolerances

Meshes.jl uses a fix absolute tolerance, which is used when comparing numbers, e.g., when using `point in geometry`. The default value
for `Float64` computations is `1e-10`, while the default for `Float32` is `1f-5`. If you want to use a custom value for the
absolute tolerance, you can use ScopedValues.jl to set them. For this you can run

```julia
using Meshes, ScopedValues
with(Meshes.ATOL64 => 1e-9, Meshes.ATOL32 => 1f-4) do
    # do your computations with the adjusted tolerances here
end
```

to set the tolerance for the `Float32` computations within in `do ... end` block to `1f-4` and for the `Float64` computations to `1e-9`.
