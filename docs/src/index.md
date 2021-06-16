# Meshes.jl

*Computational geometry and meshing algorithms in Julia.*

[![Build Status](https://img.shields.io/github/workflow/status/JuliaGeometry/Meshes.jl/CI?style=flat-square)](https://github.com/JuliaGeometry/Meshes.jl/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/JuliaGeometry/Meshes.jl?style=flat-square)](https://codecov.io/gh/JuliaGeometry/Meshes.jl)
[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue?style=flat-square)](https://JuliaGeometry.github.io/Meshes.jl/stable)
[![Latest Documentation](https://img.shields.io/badge/docs-latest-blue?style=flat-square)](https://JuliaGeometry.github.io/Meshes.jl/dev)
[![License File](https://img.shields.io/badge/license-MIT-blue?style=flat-square)](https://github.com/JuliaGeometry/Meshes.jl/blob/master/LICENSE)

## Overview

[Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl) provides efficient
implementations of concepts from computational geometry and finite element
analysis. It promotes rigorous mathematical definitions of spatial discretizations
(a.k.a. meshes) that are adequate for describing general manifolds embedded in $\R^n$,
including surfaces described with spherical coordinates, and geometries described
with multiple coordinate reference systems. Our ambitious goal is to provide all the
features of the [CGAL](https://www.cgal.org) project in pure Julia.

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

## Installation

Get the latest stable release with Julia's package manager:

```julia
] add Meshes
```

## Quick example

Although we didn't have time to document the functionality of the package properly,
we still would like to illustrate some important features. For more information on
available functionality, please consult the [Reference guide](points.md) and the
[suite of tests](https://github.com/JuliaGeometry/Meshes.jl/tree/master/test) in
the package.

```@example overview
using Meshes, MeshViz
import CairoMakie
```

### Points and vectors

A [`Point`](@ref) is defined by its coordinates in a global reference system. The type of the
coordinates is determined automatically based on the specified literals, or is forced
to a specific type using helper constructors (e.g. `Point2`, `Point3`, `Point2f`, `Point3f`).

A vector [`Vec`](@ref) follows the same pattern. It can be constructed with the generic `Vec`
constructor or with the variants `Vec2` and `Vec3` for double precision and `Vec2f`
and `Vec3f` for single precision.

```@example overview
A = Point(0, 0) # point with integer coordinates
B = Point(1, 0) # another point in 2D space
C = Point(0.0, 1.0) # double precision
D = Point2(0, 1) # double precision from Int literal
E = Point(1, 2, 3) # a point in 3D space
F = Point3(1, 2, 3) # another point now with double precision
G = Point(1f0, 2f0, 3f0) # single precision
H = Point3f(1, 2, 3) # single precision from Int literal

for P in (A,B,C,D,E,F,G,H)
  println("Coordinate type: ", coordtype(P))
  println("Embedding dimension: ", embeddim(P))
end
```

Points can be subtracted to produce a vector:

```@example overview
B - A
```

They can't be added, but their coordinates can:

```@example overview
coordinates(G) + coordinates(H)
```

We can add a point to a vector though, and get a new point:

```@example overview
G + Vec3f(1,1,1)
```

And finally, we can create points at random with:

```@example overview
ps = rand(Point2, 10)
```

### Primitives

Primitive geometries such as [`Box`](@ref), [`Ball`](@ref), [`Sphere`](@ref),
[`Cylinder`](@ref) are those geometries that can be efficiently represented
in a computer without discretization. We can construct such geometries using
clean syntax:

```@example overview
b = Box((0,0), (1,1))
```

```@example overview
s = Sphere((0,0), 1)
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
collect(vs)
```

### Polytopes

Polytopes are geometries with "flat" sides. They generalize polygons and polyhedra.
Most commonly used polytopes are already defined in the project, including
[`Segment`](@ref), [`Ngon`](@ref) (e.g. Triangle, Quadrangle), [`Tetrahedron`](@ref),
[`Pyramid`](@ref) and [`Hexahedron`](@ref).

```@example overview
t = Triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
```

Some of these geometries have additional functionality like the measure (or area):

```@example overview
measure(t) == area(t) == 1/2
```

Or the ability to know whether or not a point is inside:

```@example overview
p = Point(0.5, 0.0)

p ∈ t
```

For line segments, for example, we have robust intersection algorithms:

```@example overview
s1 = Segment((0.0,0.0), (1.0,0.0))
s2 = Segment((0.5,0.0), (2.0,0.0))

s1 ∩ s2
```

Polytopes are widely used in GIS software under names such as "LineString" and "Polygon".
We provide robust implementations of these concepts, which are formally known as polygonal
[`Chain`](@ref) and [`PolyArea`](@ref).

We can compute the orientation of a chain as clockwise or counter-clockwise, can open and
close the chain, create bridges between the various inner rings with the outer ring, and
other useful functionality:

```@example overview
p = PolyArea(Point2[(0,0), (2,0), (2,2), (1,3), (0,2), (0,0)])
```

The orientation of the above polygonal area is counter-clockwise (CCW):

```@example overview
orientation(p)
```

We can get the outer polygonal chain, and reverse it:

```@example overview
c = chains(p)[1]

reverse(c)
```

A closed chain (a.k.a. ring) has circular vertices:

```@example overview
v = vertices(c)
```

This means that we can index the vertices with indices that go
beyond the range `1:nvertices(c)`. This is very useful for
writing algorithms:

```@example overview
v[1:10]
```

We can also compute angles of any given chain, no matter if it
is open or closed:

```@example overview
angles(c) * 180 / pi
```

The sign of these angles is a function of the orientation:

```@example overview
angles(reverse(c)) * 180 / pi
```

In case of closed chains, we can compute inner angles as well:

```@example overview
innerangles(c) * 180 / pi
```

And there is a lot more functionality available like for instance
determining whether or not a polygonal area or chain is simple:

```@example overview
issimple(p)
```

### Meshes

Efficient (lazy) mesh representations are provided, including
[`CartesianGrid`](@ref) and [`SimpleMesh`](@ref):

```@example overview
g = CartesianGrid(100, 100)
```

No memory is allocated:

```@example overview
@allocated CartesianGrid(10000, 10000, 10000)
```

but we can still loop over the elements, which are quadrangles in 2D:

```@example overview
collect(elements(g))
```

We can construct a general unstructured mesh with a global vector of points
and a collection of [`Connectivity`](@ref) that store the indices to the
global vector of points:

```@example overview
points = Point2[(0,0), (1,0), (0,1), (1,1), (0.25,0.5), (0.75,0.5)]
tris  = connect.([(1,5,3),(4,6,2)], Triangle)
quads = connect.([(1,2,6,5),(4,3,5,6)], Quadrangle)
mesh = SimpleMesh(points, [tris; quads])
```

The actual geometries of the elements are materialized in a lazy fashion
like with the Cartesian grid:

```@example overview
collect(elements(mesh))
```

and all geometries and meshes can be visualized with
[MeshViz.jl](https://github.com/JuliaGeometry/MeshViz.jl):

```@example overview
viz(mesh, showfacets = true)
```

### Mesh data

To attach data to the geometries of a mesh, we can use the
[`meshdata`](@ref) function, which combines a mesh object
with a collection of Tables.jl tables. For example, it is
common to attach a table `vtable` to the vertices and a
table `etable` to the elements of the mesh:

```@example overview
d = meshdata(mesh,
  vtable = (temperature=rand(6), pressure=rand(6)),
  etable = (quality=["A","B"], state=[true,false])
)
```

More generally, we can attach a table to any rank:

- 0 (vertices)
- 1 (segments)
- 2 (triangles, quadrangles, ...)
- 3 (tetrahedrons, hexahedrons, ...)

To retrieve the data table for a given rank we use
the `values` function:

```@example overview
values(d, 0)
```

```@example overview
values(d, 2)
```

If we ommit the rank, the function will return the `etable`
of the mesh:

```@example overview
values(d)
```

When a table is not available for a given rank, the value
`nothing` is returned instead:

```@example overview
values(d, 1) === nothing
```

Finally, we can use the `domain` function to retrieve the
underlying domain of the data, which in this case is a
`SimpleMesh`:

```@example overview
domain(d)
```
