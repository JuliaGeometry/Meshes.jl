# Meshes.jl

*Computational geometry and meshing algorithms in Julia.*

[![Build Status](https://img.shields.io/github/actions/workflow/status/JuliaGeometry/Meshes.jl/CI.yml?branch=master&style=flat-square)](https://github.com/JuliaGeometry/Meshes.jl/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/JuliaGeometry/Meshes.jl?style=flat-square)](https://codecov.io/gh/JuliaGeometry/Meshes.jl)
[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue?style=flat-square)](https://JuliaGeometry.github.io/Meshes.jl/stable)
[![Latest Documentation](https://img.shields.io/badge/docs-latest-blue?style=flat-square)](https://JuliaGeometry.github.io/Meshes.jl/dev)
[![License File](https://img.shields.io/badge/license-MIT-blue?style=flat-square)](https://github.com/JuliaGeometry/Meshes.jl/blob/master/LICENSE)

## Overview

[Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl) provides efficient
implementations of concepts from computational geometry. It promotes rigorous
mathematical definitions of spatial discretizations (a.k.a. meshes) that are
adequate for describing general manifolds embedded in $\R^n$, including surfaces
described with spherical coordinates, and geometries described with multiple
coordinate reference systems.

Unlike other existing efforts in the Julia ecosystem, this project is being carefully
designed to facilitate the use of *meshes across different scientific domains*. We
follow a strict set of good software engineering practices, and are quite pedantic
in our test suite to make sure that all our implementations are free of bugs in both
single and double floating point precision. Additionally, we guarantee type stability.

The design of this project was motivated by various issues encountered with past attempts
to represent geometry, which have been originally designed for visualization
purposes (e.g. [GeometryTypes.jl](https://github.com/JuliaGeometry/GeometryTypes.jl),
[GeometryBasics.jl](https://github.com/JuliaGeometry/GeometryBasics.jl)) or specifically
for finite element analysis (e.g. [JuAFEM.jl](https://kristofferc.github.io/JuAFEM.jl/dev/manual/grid),
[MeshCore.jl](https://github.com/PetrKryslUCSD/MeshCore.jl)). We hope to provide a smoother
experience with mesh representations that are adequate for finite finite element analysis,
advanced geospatial modeling *and* visualization, not just one domain.

For advanced data science with geospatial data (i.e., tables over meshes), consider the
[GeoStats.jl](https://github.com/JuliaEarth/GeoStats.jl) framework. It provides sophisticated
methods for estimating (interpolating), simulating and learning geospatial functions over
[Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl) meshes. Please check the
[Geospatial Data Science with Julia](https://juliaearth.github.io/geospatial-data-science-with-julia)
book for more information:


```@raw html
<p align="center">
  <a href="https://juliaearth.github.io/geospatial-data-science-with-julia">
    <img src="https://juliaearth.github.io/geospatial-data-science-with-julia/images/cover.svg" width="200px" hspace="20">
  </a>
</p>
```

If you have questions or would like to brainstorm ideas in general, don't hesitate to start
a thread in our [zulip channel](https://julialang.zulipchat.com/#narrow/stream/275558-meshes.2Ejl).
We are happy to improve the ecosystem to meet user's needs.

## Installation

Get the latest stable release with Julia's package manager:

```
] add Meshes
```

## Quick example

Although we didn't have time to document the functionality of the package properly,
we still would like to illustrate some important features. For more information on
available functionality, please consult the [Reference guide](vectors.md) and the
[suite of tests](https://github.com/JuliaGeometry/Meshes.jl/tree/master/test) in
the package.

In all examples we assume the following packages are loaded:

```@example overview
using Meshes
import CairoMakie as Mke
```

### Points and vectors

A [`Point`](@ref) is defined by its coordinates in a global reference system. The type of the
coordinates is determined automatically based on the specified literals.
`Integer` coordinates are converted to `Float64` to fulfill the requirements of most
geometric processing algorithms, which would be undefined in a discrete scale.

A vector [`Vec`](@ref) follows the same pattern. It can be constructed with the `Vec` constructor.

```@example overview
Point(0.0, 1.0) # double precision as expected
```

```@example overview
Point(0f0, 1f0) # single precision as expected
```

```@example overview
Point(0, 0) # Integer is converted to Float64 by design
```

```@example overview
Point(1.0, 2.0, 3.0) # double precision as expected
```

```@example overview
Point(1f0, 2f0, 3f0) # single precision as expected
```

```@example overview
Point(1, 2, 3) # Integer is converted to Float64 by design
```

Points can be subtracted to produce a vector:

```@example overview
A = Point(1.0, 1.0)
B = Point(3.0, 3.0)
B - A
```

They can't be added, but their coordinates can:

```@example overview
coordinates(A) + coordinates(B)
```

We can add a point to a vector though, and get a new point:

```@example overview
A + Vec(1, 1)
```

And finally, we can create points at random with:

```@example overview
rand(Point{2})
```

### Primitives

Primitive geometries such as [`Box`](@ref), [`Ball`](@ref), [`Sphere`](@ref),
[`Cylinder`](@ref) are those geometries that can be efficiently represented
in a computer without discretization. We can construct such geometries using
clean syntax:

```@example overview
b = Box((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))

viz(b)
```

```@example overview
s = Sphere((0.0, 0.0, 0.0), 1.0)

viz(s)
```

The parameters of these primitive geometries can be queried easily:

```@example overview
extrema(b)
```

```@example overview
centroid(s), radius(s)
```

As well as their measure (e.g. area, volume) and other geometric properties:

```@example overview
measure(b)
```

We can sample random points on primitives using different methods:

```@example overview
vs = sample(s, RegularSampling(10)) # 10 points over the sphere
```

And collect the generator with:

```@example overview
viz(collect(vs))
```

### Polytopes

Polytopes are geometries with "flat" sides. They generalize polygons and polyhedra.
Most commonly used polytopes are already defined in the project, including
[`Segment`](@ref), [`Ngon`](@ref) (e.g. Triangle, Quadrangle), [`Tetrahedron`](@ref),
[`Pyramid`](@ref) and [`Hexahedron`](@ref).

```@example overview
t = Triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))

viz(t)
```

Some of these geometries have additional functionality like the measure (or area):

```@example overview
measure(t)
```

```@example overview
measure(t) == area(t)
```

Or the ability to know whether or not a point is inside:

```@example overview
p = Point(0.5, 0.0)

p ∈ t
```

For line segments, for example, we have robust intersection algorithms:

```@example overview
s1 = Segment((0.0, 0.0), (1.0, 0.0))
s2 = Segment((0.5, 0.0), (2.0, 0.0))

s1 ∩ s2
```

Polytopes are widely used in GIS software under names such as "LineString" and "Polygon".
We provide robust implementations of these concepts, which are formally known as polygonal
[`Chain`](@ref) and [`PolyArea`](@ref).

We can compute the orientation of a chain as clockwise or counter-clockwise, can open and
close the chain, create bridges between the various inner rings with the outer ring, and
other useful functionality:

```@example overview
p = PolyArea((0,0), (2,0), (2,2), (1,3), (0,2))

viz(p)
```

The orientation of the above polygonal area is counter-clockwise (CCW):

```@example overview
orientation(p)
```

We can get the outer ring, and reverse it:

```@example overview
r = rings(p) |> first

reverse(r)
```

A ring has circular vertices:

```@example overview
v = vertices(r)
```

This means that we can index the vertices with indices that go
beyond the range `1:nvertices(r)`. This is very useful for
writing algorithms:

```@example overview
v[1:10]
```

We can also compute angles of any given chain, no matter if it
is open or closed:

```@example overview
angles(r) * 180 / pi
```

The sign of these angles is a function of the orientation:

```@example overview
angles(reverse(r)) * 180 / pi
```

In case of rings (i.e. closed chains), we can compute inner angles as well:

```@example overview
innerangles(r) * 180 / pi
```

And there is a lot more functionality available like for instance
determining whether or not a polygonal area or chain is simple:

```@example overview
issimple(p)
```

### Meshes

Efficient (lazy) mesh representations are provided, including
[`CartesianGrid`](@ref) and [`SimpleMesh`](@ref), which are
specific types of [`Domain`](@ref):

```@docs
Domain
```

```@example overview
grid = CartesianGrid(100, 100)

viz(grid, showsegments = true)
```

No memory is allocated:

```@example overview
@allocated CartesianGrid(10000, 10000, 10000)
```

but we can still loop over the elements, which are quadrangles in 2D:

```@example overview
collect(grid)
```

We can construct a general unstructured mesh with a global vector of points
and a collection of [`Connectivity`](@ref) that store the indices to the
global vector of points:

```@example overview
points = [(0,0), (1,0), (0,1), (1,1), (0.25,0.5), (0.75,0.5)]
tris  = connect.([(1,5,3), (4,6,2)], Triangle)
quads = connect.([(1,2,6,5), (4,3,5,6)], Quadrangle)
mesh = SimpleMesh(points, [tris; quads])
```

```@example overview
viz(mesh, showsegments = true)
```

The actual geometries of the elements are materialized in a lazy fashion
like with the Cartesian grid:

```@example overview
collect(mesh)
```