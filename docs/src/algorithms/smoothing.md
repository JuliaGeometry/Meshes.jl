# Smoothing

```@docs
smooth
SmoothingMethod
```

## Taubin

```@docs
TaubinSmoothing
```

```@example taubin
using Meshes
using MeshViz
using PlyIO

import WGLMakie as Mke

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
smesh = smooth(mesh, TaubinSmoothing(30))

fig = Mke.Figure(resolution = (800, 1200))
viz(fig[1,1], mesh)
viz(fig[2,1], smesh)
fig
```
