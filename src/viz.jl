# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    viz(object; [options])

Visualize Meshes.jl `object` with various `options`.

## Available options

* `color`        - color of geometries
* `alpha`        - transparency in [0,1]
* `colormap`     - color scheme/map from ColorSchemes.jl
* `colorrange`   - minimum and maximum color values
* `showsegments` - visualize segments
* `segmentcolor` - color of segments
* `segmentsize`  - width of segments
* `showpoints`   - visualize points
* `pointcolor`   - color of points
* `pointsize`    - size of points

The option `color` can be a single scalar or a vector
of scalars. For [`Mesh`](@ref) subtypes, the length of
the vector of colors determines if the colors should be
assigned to vertices or to elements.

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

* This function will only work in the presence of
  a Makie.jl backend via package extensions in
  Julia v1.9 or later versions of the language.
"""
function viz end

"""
    viz!(object; [options])

Visualize Meshes.jl `object` in an existing
scene with `options` forwarded to [`viz`](@ref).
"""
function viz! end
