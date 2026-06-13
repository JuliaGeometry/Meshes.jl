# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    viz(object; [options])

Visualize Meshes.jl `object` (e.g., mesh, geometry)
with `options` such as `color` and `alpha` values for
transparency. All available options will be documented
below upon loading a Makie.jl backend.

## Examples

Different coloring methods (vertex vs. element):

```
# vertex coloring (i.e. linear interpolation)
viz(mesh, color = 1:nvertices(mesh))

# element coloring (i.e. discrete colors)
viz(mesh, color = 1:nelements(mesh))
```

Different strategies to show the boundary of
geometries (showsegments vs. boundary):

```
# visualize boundary with showsegments
viz(polygon, showsegments = true)

# visualize boundary with separate call
viz(polygon)
viz!(boundary(polygon))
```

### Notes

This function will only work in the presence of a Makie.jl
backend via package extensions in Julia v1.9 or later
versions of the language.

Values passed to the `color` option can be of any type, as
long as they can be converted to colors by the `colorfy`
function from the Colorfy.jl package.

For [`Mesh`](@ref) subtypes, the length of the vector of
colors determines if the colors should be assigned to
vertices or to elements.
"""
function viz end

"""
    viz!(object; [options])

Visualize Meshes.jl `object` (e.g., mesh, geometry) in an
existing scene with `options` forwarded to [`viz`](@ref).
"""
function viz! end
