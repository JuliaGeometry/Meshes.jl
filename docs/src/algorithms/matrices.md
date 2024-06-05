# Matrices

```@example matrices
using Meshes # hide
import CairoMakie as Mke # hide
```

## Laplace

```@docs
laplacematrix
```

```@example matrices
grid = CartesianGrid(10, 10)

laplacematrix(grid, kind = :uniform)
```

```@example matrices
points = [(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)]
connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)])
mesh = SimpleMesh(points, connec)

laplacematrix(mesh, kind = :cotangent)
```

## Measure

```@docs
measurematrix
```

```@example matrices
grid = CartesianGrid(10, 10)

measurematrix(grid)
```

```@example matrices
points = [(0, 0), (1, 0), (0, 1), (1, 1), (0.5, 0.5)]
connec = connect.([(1, 2, 5), (2, 4, 5), (4, 3, 5), (3, 1, 5)])
mesh = SimpleMesh(points, connec)

measurematrix(mesh)
```

## Adjacency

```@docs
adjacencymatrix
```

```@example matrices
grid = CartesianGrid(10, 10)

adjacencymatrix(grid)
```

```@example matrices
adjacencymatrix(grid, rank = 0)
```
