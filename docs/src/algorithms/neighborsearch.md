# Neighbor search

```@example neighbors
using Meshes # hide
import CairoMakie as Mke # hide
```

It is often useful to search neighbor elements in a domain given a
point of reference. This can be performed with search methods:

```@docs
NeighborSearchMethod
BoundedNeighborSearchMethod
search
search!
searchdists
searchdists!
```

Search methods are constructed with various types of parameters.
One may be interested in k-nearest neighbors, or interested in
neighbors within a certain [`Neighborhood`](@ref):

```@docs
Neighborhood
MetricBall
```

The following example demonstrates neighbor search with the
[`KNearestSearch`](@ref) method:

```@example neighbors
grid = CartesianGrid(10, 10)

# 4-nearest neighbors
searcher = KNearestSearch(grid, 4)

inds = search(Point(5.0, 5.0), searcher)
```

The function [`search`](@ref) returns the indices of the elements
in the domain that are neighbors of the point. The elements are:

```@example neighbors
grid[inds]
```

Alternatively, the function [`searchdists`](@ref) also returns
the distances to the (centroids) of the elements:

```@example neighbors
inds, dists = searchdists(Point(5.0, 5.0), searcher)

dists
```

Finally, the functions [`search!`](@ref) and [`searchdists!`](@ref)
can be used in hot loops to avoid unnecessary memory allocations.

## BallSearch

```@docs
BallSearch
```

## KNearestSearch

```@docs
KNearestSearch
```

## KBallSearch

```@docs
KBallSearch
```
