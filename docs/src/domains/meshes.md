# Meshes

```@example meshes
using Meshes # hide
import CairoMakie as Mke # hide
```

Meshes can be constructed directly (e.g. [`CartesianGrid`](@ref)) or based on other
constructs such as connectivity lists and topological structures (e.g. [`SimpleMesh`](@ref)).

## Overview

```@docs
Mesh
Grid
```

```@docs
CartesianGrid
```

```@example meshes
# 3D Cartesian grid
grid = CartesianGrid(10, 10, 10)

viz(grid, showsegments = true)
```

```@docs
RectilinearGrid
```

```@example meshes
# 2D rectilinear grid
x = 0.0:0.2:1.0
y = [0.0, 0.1, 0.3, 0.7, 0.9, 1.0]
grid = RectilinearGrid(x, y)

viz(grid, showsegments = true)
```

```@docs
StructuredGrid
```

```@example meshes
# 2D structured grid
X = [i/20 * cos(3π/2 * (j-1) / (30-1)) for i in 1:20, j in 1:30]
Y = [i/20 * sin(3π/2 * (j-1) / (30-1)) for i in 1:20, j in 1:30]
grid = StructuredGrid(X, Y)

viz(grid, showsegments = true)
```

```@docs
SimpleMesh
```

```@example meshes
# global vector of 2D points
points = [(0,0),(1,0),(0,1),(1,1),(0.25,0.5),(0.75,0.5)]

# connect the points into N-gon
connec = connect.([(1,2,6,5),(2,4,6),(4,3,5,6),(3,1,5)], Ngon)

# 2D mesh made of N-gon elements
mesh = SimpleMesh(points, connec)

viz(mesh, showsegments = true)
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
GridTopology
HalfEdgeTopology
SimpleTopology
```

## Relations

```@docs
TopologicalRelation
Boundary
Coboundary
Adjacency
```

Consider the following examples with the [`Boundary`](@ref) and
[`Coboundary`](@ref) relations defined for the [`HalfEdgeTopology`](@ref):

```@example meshes
# global vector of 2D points
points = [(0,0),(1,0),(0,1),(1,1),(0.25,0.5),(0.75,0.5)]

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

## Matrices

Based on topological relations, we can extract matrices that
are widely used in applications such as [`laplacematrix`](@ref),
and [`adjacencymatrix`](@ref).

### Laplace

```@docs
laplacematrix
```

```@example meshes
grid = CartesianGrid(10, 10)

laplacematrix(grid, kind = :uniform)
```

```@example meshes
points = [(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)]
connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)])
mesh = SimpleMesh(points, connec)

laplacematrix(mesh, kind = :cotangent)
```

### Measure

```@docs
measurematrix
```

```@example meshes
grid = CartesianGrid(10, 10)

measurematrix(grid)
```

```@example meshes
points = [(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)]
connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)])
mesh = SimpleMesh(points, connec)

measurematrix(mesh)
```

### Adjacency

```@docs
adjacencymatrix
```

```@example meshes
grid = CartesianGrid(10, 10)

adjacencymatrix(grid)
```

```@example meshes
adjacencymatrix(grid, rank = 0)
```
