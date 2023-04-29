# Transforms

```@example transforms
using JSServe: Page # hide
Page(exportable=true, offline=true) # hide
```

```@example transforms
using Meshes, MeshViz # hide
import WGLMakie as Mke # hide
```

Geometric (e.g. coordinates) transforms are implemented according to the
[TransformsBase.jl](https://github.com/JuliaML/TransformsBase.jl) interface.
Please read their documentation for more details.

```@docs
GeometricTransform
StatelessGeometricTransform
```

## Rotate

```@docs
Rotate
```

```@example transforms
grid = CartesianGrid(10, 10)

mesh = grid |> Rotate(ClockwiseAngle(π/4))

fig = Mke.Figure(resolution = (800, 400))
viz(fig[1,1], grid)
viz(fig[1,2], mesh)
fig
```

## Translate

```@docs
Translate
```

```@example transforms
grid = CartesianGrid(10, 10)

mesh = grid |> Translate(10., 20.)

fig = Mke.Figure(resolution = (800, 400))
viz(fig[1,1], grid)
viz(fig[1,2], mesh)
fig
```

## Stretch

```@docs
Stretch
```

```@example transforms
grid = CartesianGrid(10, 10)

mesh = grid |> Stretch(2., 3.)

fig = Mke.Figure(resolution = (800, 400))
viz(fig[1,1], grid)
viz(fig[1,2], mesh)
fig
```

## StdCoords

```@docs
StdCoords
```

```@example transforms
# Cartesian grid with coordinates [0,10] x [0,10]
grid = CartesianGrid(10, 10)

# scale coordinates to [-1,1] x [-1,1]
mesh = grid |> StdCoords()

fig = Mke.Figure(resolution = (800, 400))
viz(fig[1,1], grid)
viz(fig[1,2], mesh)
fig
```

## Repair

```@docs
Repair
```

```@example transforms
# mesh with unreferenced point
points = Point3[(0, 0, 0), (0, 0, 1), (5, 5, 5), (0, 1, 0), (1, 0, 0)]
connec = connect.([(1, 2, 4), (1, 2, 5), (1, 4, 5), (2, 4, 5)])
mesh   = SimpleMesh(points, connec)

rmesh = mesh |> Repair{1}()
```

## Smoothing

```@docs
LambdaMuSmoothing
LaplaceSmoothing
TaubinSmoothing
```

```@example transforms
using PlyIO

# helper function to read *.ply files
function readply(fname)
  ply = load_ply(fname)
  x = ply["vertex"]["x"]
  y = ply["vertex"]["y"]
  z = ply["vertex"]["z"]
  points = Point3.(x, y, z)
  connec = [connect(Tuple(c.+1)) for c in ply["face"]["vertex_indices"]]
  SimpleMesh(points, connec)
end

# download mesh from the web
file = download(
  "https://raw.githubusercontent.com/juliohm/JuliaCon2021/master/data/beethoven.ply"
)

# read mesh from disk
mesh = readply(file)

# smooth mesh with 30 iterations
smesh = mesh |> TaubinSmoothing(30)

fig = Mke.Figure(resolution = (800, 1200))
viz(fig[1,1], mesh)
viz(fig[2,1], smesh)
fig
```
