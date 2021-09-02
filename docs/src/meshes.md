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

```@example meshes
# convert topology to half-edge topology for fast lookup
ht = convert(HalfEdgeTopology, topology(mesh))
# We can create the boundary facets of the first Ngon.
# The boundary is parametrized on the dimension we are mapping from,
# and the lower dimension we are mapping to.
# The Ngons are of order 2 and are made up of segments with order 1:
b = Boundary{2,1}(ht)
# Show boundary of Ngon number one:
b(1)
```

```@example meshes
## The co-boundary provides the connectivity from a lower to higher
# dimension.
cb = Coboundary{1, 2}(ht)
# From this, we can look up the pair of Ngons that share a given 
# segment in their boundary.
cb(3)
```

```@example meshes
# For segments that are on the boundary of the domain itself, we
# only see a single segment:
cb(1)
```

```@example meshes
## We can loop over the facets and extract geometrical properties:
for (i, f) in enumerate(facets(ht))
    seg = materialize(f, points)
    n = cb(i)
    if length(n) == 2
        println("$i: Interior facet:\n - connected to Ngons: $(n)")
    else
        println("$i: Exterior facet:\n - connected to Ngon: $(n[1])")
    end
    l = measure(seg)
    println(" - segment length = $l")
    c = centroid(seg)
    println(" - centroid at $c")
end
```
