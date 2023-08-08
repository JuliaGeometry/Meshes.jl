# Trajectories

```@example trajec
using JSServe: Page # hide
Page(exportable=true, offline=true) # hide
```

```@example trajec
using Meshes # hide
import WGLMakie as Mke # hide
```

Trajectories of geometries are special geometry sets
with one-dimensional grid topology. They are often
used in geosciences to represent drill holes, wells,
etc.

```@docs
CylindricalTrajectory
```

```@example trajec
# construct centroids along Bezier curve
b = BezierCurve([(0, 0, 0), (3, 3, 0), (3, 0, 7)])
c = [b(t) for t in range(0, stop=1, length=20)]

# cylindrical trajectory with radius 2
CylindricalTrajectory(c, 2) |> viz
```
