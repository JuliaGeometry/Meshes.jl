# Preferences

## Setting the absolute tolerance

Meshes.jl uses a fix absolute tolerance, which is used when comparing numbers, e.g., when using `point in geometry`. The default value
for `Float64` computations is `1.0e-10`, while the default for `Float32` is `1.0f-5`. If you want to use a custom value for the
absolute tolerance, you can use Preferences.jl to set them. For this you can run

```julia
using Meshes, Preferences
set_preferences!(Meshes, "atol_float32" => 1.0f-4)
```

to set the tolerance for the `Float32` computations to `1.0f-4`. Analogously, you can set the preference `atol_float64`.
Calling `set_preferences!` will create a file "LocalPreferences.toml", which stores the preferences. After changing a
preference, you will need to restart the REPL for the changes to take effect.
To switch back to the default values, you can simply delete the "LocalPreferences.toml" file or call `delete_preferences!(Meshes, "atol_float32")`.
