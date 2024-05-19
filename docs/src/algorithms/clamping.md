# Clamping

Meshes adds methods to Julia's built-in `clamp` function. The additional methods clamp points to the edges of a box in any number of dimensions. The target points and boxes must have the same number of dimensions and the same numeric type.

```@docs
clamp(::Point{Dim}, ::Box{Dim}) where {Dim}
clamp(::PointSet, ::Box)
```

```@example clamping
using Meshes # hide
import CairoMakie as Mke # hide

# set of 2D points to clamp
points = PointSet(rand(2, 100))

# 2D box defining the clamping boundaries
box = Box((0.25, 0.25), (0.75, 0.75))

# clamp point coordinates to the box edges
clamped = clamp(points, box)

fig = Mke.Figure(size=(800, 400))
ax = Mke.Axis(fig[1,1], title="unclamped", aspect=1, limits=(0,1,0,1))
viz!(ax, box)
viz!(ax, points, color=:black, pointsize=6)
ax = Mke.Axis(fig[1,2], title="clamped", aspect=1, limits=(0,1,0,1))
viz!(ax, box)
viz!(ax, clamped, color=:black, pointsize=6)
fig
```
