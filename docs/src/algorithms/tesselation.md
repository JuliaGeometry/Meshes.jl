# Tesselation

```@example tesselation
using Meshes # hide
import CairoMakie as Mke # hide
```

```@docs
tesselate
TesselationMethod
```

## DelaunayTesselation

```@docs
DelaunayTesselation
```

```@example tesselation
points = rand(Point{2}, 100)

mesh = tesselate(points, DelaunayTesselation())

viz(mesh, showsegments = true)
viz!(points, color = :red)
Mke.current_figure()
```

## VoronoiTesselation

```@docs
VoronoiTesselation
```

```@example tesselation
points = rand(Point{2}, 100)

mesh = tesselate(points, VoronoiTesselation())

viz(mesh, showsegments = true)
viz!(points, color = :red)
Mke.current_figure()
```
