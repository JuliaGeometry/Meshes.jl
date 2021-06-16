# Meshes

## Overview

Meshes can be constructed directly (e.g. [`CartesianGrid`](@ref)) or based on other
constructs such as connectivity lists and topological structures (e.g. [`SimpleMesh`](@ref)).

```@docs
Mesh
CartesianGrid
SimpleMesh
```

### Connectivities

```@docs
Connectivity
connect
materialize
```

### Topology

```@docs
Topology
FullTopology
GridTopology
HalfEdgeTopology
```

## Examples

```@example meshes
using Meshes

# 3D Cartesian grid
grid = CartesianGrid(10, 10, 10)
```

```@example meshes
using MeshViz
import CairoMakie

viz(grid, showfacets = true)
```

```@example meshes
# global vector of 2D points
points = Point2[(0,0),(1,0),(0,1),(1,1),(0.25,0.5),(0.75,0.5)]

# connect the points into N-gon
connec = connect.([(1,2,6,5),(2,4,6),(4,3,5,6),(3,1,5)], Ngon)
```

```@example meshes
# 2D mesh made of N-gon elements
mesh = SimpleMesh(points, connec)
```

```@example meshes
viz(mesh, showfacets = true)
```
