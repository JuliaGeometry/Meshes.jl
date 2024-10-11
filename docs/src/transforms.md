# Transforms

```@example transforms
using Meshes # hide
using CoordRefSystems # hide
import CairoMakie as Mke # hide
```

Geometric (e.g. coordinates) transforms are implemented according to the
[TransformsBase.jl](https://github.com/JuliaML/TransformsBase.jl) interface.
Please read their documentation for more details.

```@docs
GeometricTransform
CoordinateTransform
```

Some transforms have an inverse that can be created with the [`inverse`](@ref) function.
The function [`isinvertible`](@ref) can be used to check if a transform is invertible.

```@docs
inverse
isinvertible
```

## Rotate

```@docs
Rotate
```

```@example transforms
grid = CartesianGrid(10, 10)

mesh = grid |> Rotate(π/4)

fig = Mke.Figure(size = (800, 400))
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

fig = Mke.Figure(size = (800, 400))
viz(fig[1,1], grid)
viz(fig[1,2], mesh)
fig
```

## Scale

```@docs
Scale
```

```@example transforms
grid = CartesianGrid(10, 10)

mesh = grid |> Scale(2., 3.)

fig = Mke.Figure(size = (800, 400))
viz(fig[1,1], grid)
viz(fig[1,2], mesh)
fig
```

## Affine

```@docs
Affine
```

```@example transforms
using Rotations: Angle2d

grid = CartesianGrid(10, 10)

mesh = grid |> Affine(Angle2d(π/4), [10., 20.])

fig = Mke.Figure(size = (800, 400))
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

fig = Mke.Figure(size = (800, 400))
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

fig = Mke.Figure(size = (800, 400))
viz(fig[1,1], grid)
viz(fig[1,2], mesh)
fig
```

## Proj

```@docs
Proj
```

```@example transforms
# load coordinate reference system
using CoordRefSystems: Polar

# triangle with Cartesian coordinates
triangle = Triangle((0, 0), (1, 0), (1, 1))

# reproject to polar coordinates
triangle |> Proj(Polar)
```

## Morphological

```@docs
Morphological
```

```@example transforms
# triangle with Cartesian coordinates
triangle = Triangle((0, 0), (1, 0), (1, 1))

# transform triangle coordinates
triangle |> Morphological(c -> Cartesian(c.x, c.y, zero(c.x)))
```

## LengthUnit

```@docs
LengthUnit
```

```@example transforms
using Unitful: m, cm

# convert meters to centimeters
Point(1m, 2m, 3m) |> LengthUnit(cm)
```

## Shadow

```@docs
Shadow
```

```@example transforms
ball = Ball((0, 0, 0), 1)
disk = ball |> Shadow("xy")


fig = Mke.Figure(size = (800, 400))
viz(fig[1,1], ball)
viz(fig[1,2], disk)
fig
```

## Slice

```@docs
Slice
```

```@example transforms
grid = CartesianGrid(10, 10)
subgrid = grid |> Slice(x=(1.5, 6.5), y=(3.5, 8.5))

fig = Mke.Figure(size = (800, 400))
viz(fig[1,1], grid)
viz(fig[1,2], subgrid)
fig
```

## Repair

```@docs
Repair
```

```@example transforms
# mesh with unreferenced point
points = [(0, 0, 0), (0, 0, 1), (5, 5, 5), (0, 1, 0), (1, 0, 0)]
connec = connect.([(1, 2, 4), (1, 2, 5), (1, 4, 5), (2, 4, 5)])
mesh   = SimpleMesh(points, connec)

rmesh = mesh |> Repair(1)
```

## Bridge

```@docs
Bridge
```

```@example transforms
# polygon with two holes
outer = [(0, 0), (1, 0), (1, 1), (0, 1)]
hole1 = [(0.2, 0.2), (0.2, 0.4), (0.4, 0.4), (0.4, 0.2)]
hole2 = [(0.6, 0.2), (0.6, 0.4), (0.8, 0.4), (0.8, 0.2)]
poly = PolyArea([outer, hole1, hole2])

# polygon with single outer ring
bpoly = poly |> Bridge(0.01)

fig = Mke.Figure(size = (800, 400))
viz(fig[1,1], poly)
viz(fig[1,2], bpoly)
fig
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
  points = Point.(x, y, z)
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

fig = Mke.Figure(size = (800, 1200))
viz(fig[1,1], mesh)
viz(fig[2,1], smesh)
fig
```
