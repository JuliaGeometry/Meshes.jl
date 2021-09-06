# Meshes

Meshes can be constructed directly (e.g. [`CartesianGrid`](@ref)) or based on other
constructs such as connectivity lists and topological structures (e.g. [`SimpleMesh`](@ref)).

## Overview

```@docs
Mesh
CartesianGrid
SimpleMesh
```

### Examples

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

## Connectivities

```@docs
Connectivity
connect
materialize
```

## Topology

```@docs
Topology
FullTopology
GridTopology
HalfEdgeTopology
```

## Relations

```@docs
TopologicalRelation
Boundary
Coboundary
Adjacency
```

### Examples

```@example meshes
# global vector of 2D points
points = Point2[(0,0),(1,0),(0,1),(1,1),(0.25,0.5),(0.75,0.5)]

# connect the points into N-gon
connec = connect.([(1,2,6,5),(2,4,6),(4,3,5,6),(3,1,5)], Ngon)

# 2D mesh made of N-gon elements
mesh = SimpleMesh(points, connec)
```

```@example meshes
# convert topology to half-edge topology
topo = convert(HalfEdgeTopology, topology(mesh))

# boundary relation from faces (dim=2) to edges (dim=1)
∂₂₁ = Boundary{2,1}(topo)

# show boundary of first n-gon
∂₂₁(1)
```

```@example meshes
# co-boundary relation from edges (dim=1) to faces(dim=2)
𝒞₁₂ = Coboundary{1,2}(topo)

# show n-gons that share edge 3
𝒞₁₂(3)
```