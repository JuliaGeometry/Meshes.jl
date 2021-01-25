# Meshes.jl

*Computational geometry and meshing algorithms in Julia.*

[![Build Status](https://img.shields.io/github/workflow/status/JuliaGeometry/Meshes.jl/CI?style=flat-square)](https://github.com/JuliaGeometry/Meshes.jl/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/JuliaGeometry/Meshes.jl?style=flat-square)](https://codecov.io/gh/JuliaGeometry/Meshes.jl)
[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue?style=flat-square)](https://JuliaGeometry.github.io/Meshes.jl/stable)
[![Latest Documentation](https://img.shields.io/badge/docs-latest-blue?style=flat-square)](https://JuliaGeometry.github.io/Meshes.jl/latest)
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
single and double floating point precision.

The design of this project was motivated by various conceptual issues and bugs
encountered with past attemps to represent geometry, which have been originally
designed for visualization purposes (e.g.
[GeometryTypes.jl](https://github.com/JuliaGeometry/GeometryTypes.jl),
[GeometryBasics.jl](https://github.com/JuliaGeometry/GeometryBasics.jl)).
We hope to provide a smoother experience to end users, as well as mesh representations
that are adequate for finite element analysis and advanced geospatial modeling.

## Installation

Get the latest stable release with Julia's package manager:

```julia
] add Meshes
```

## Quick example

While we didn't have time to document the functionality of the package properly,
we still would like to illustrate some of the concepts below. For more information,
we kindly ask users to read the test suite.

```@example overview
using Meshes
```

### Points and vectors

A `Point` is defined by its coordinates in a global reference system. The type of the
coordinates is determined automatically based on the specified literals, or is forced
to a specific type using helper constructors (e.g. `Point2`, `Point3`, `Point2f`, `Point3f`).

A vector `Vec` follows the same pattern. It can be constructed with the generic `Vec`
constructor or with the variants `Vec2` and `Vec3` for double precision and `Vec2f`
and `Vec3f` for single precision.

```@example overview
A = Point(0, 0) # point with integer coordinates
B = Point(1, 0) # another point in 2D space

# points can be subtracted to produce a vector
v = B - A
```

```@example overview
C = Point(0.0, 1.0) # double precision
D = Point2(0, 1) # double precision from Int literal

C == D
```

```@example overview
# 3D points follow the same pattern
E = Point(1, 2, 3)
F = Point3(1, 2, 3)

coordtype(E), coordtype(F)
```

```@example overview
# for single precision, we add a `f` to the constructors
G = Point(1f0, 2f0, 3f0) # single precision
H = Point3f(1, 2, 3) # single precision from Int literal
```

```@example overview
# points can't be added, only their coordinates
coordinates(G) + coordinates(H)
```

```@example overview
# alternatively, we can add a point to a vector
G + Vec3f(1,1,1)
```

```@example overview
# any point lives in an embedding space of given dimension
embeddim(G)
```

```@example overview
# we can create points at random
ps = rand(Point2, 10)
```

### Primitives

Primitive geometries such as `Box`, `Sphere`, `Cylinder` are those geometries
that can be efficiently represented in a computer without discretization.
We can construct such geometries using clean syntax:

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
center(s), radius(s)
```

As well as their measure (e.g. area, volume) and other geometric properties:

```@example overview
measure(b) == 1
```

### Polytopes

Polytopes are geometries with "flat" sides. They generalize polygons and polyhedra.
Most commonly used polytopes are already defined in the project, including `Triangle`
`Pyramid`, `Quadrangle`, `Segment`, `Tetrahedron`, and `Hexahedron`. Their vertices
follow the GMSH ordering convention by default, but this is also customizable.

```@example overview
t = Triangle((0,0), (1,0), (0,1))
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
`Chain` and `PolyArea`.

We can compute the orientation of a chain as clockwise or counter-clockwise, can open and
close the chain, create bridges between the various inner rings with the outer ring, and
other useful functionality:

```@example overview
p = PolyArea(Point2[(0,0), (2,0), (2,2), (1,3), (0,2), (0,0)])
```

```@example overview
# orientation is counter-clockwise (CCW)
orientation(p)
```

```@example overview
# get outer polygonal chain
c = chains(p)[1]

# revert its vertices
reverse(c)
```

```@example overview
# a closed chain (a.k.a ring) has circular vertices
v = vertices(c)
```

```@example overview
# we can index beyond 1:nvertices(c)
v[1:10]
```

```@example overview
# let's compute the angles of the chain
angles(c) * 180 / pi
```

```@example overview
# angles are a function of the orientation
angles(reverse(c)) * 180 / pi
```

```@example overview
# in the case of closed chains, we can ask for the inner angles
innerangles(c) * 180 / pi
```

And there is much more functionality available like for instance
determining whether or not a polygonal area or chain is simple:

```@example overview
issimple(p)
```

### Meshes

Efficient (lazy) mesh representations are provided to avoid unncessary
memory allocations. For example, we can create a `CartesianGrid` or
general `UnstructuredMesh`:

```@example overview
g = CartesianGrid(100, 100)
```

```@example overview
# no memory is allocated
@allocated CartesianGrid(10000, 10000, 10000)
```

```@example overview
# we can still loop over the elements
collect(elements(g))
```

```@example overview
# here is a general mesh with multiple element types
points = Point2[(0,0), (1,0), (0,1), (1,1), (0.25,0.5), (0.75,0.5)]
Δs = connect.([(3,1,5),(4,6,2)], Triangle)
□s = connect.([(1,2,5,6),(5,6,3,4)], Quadrangle)
mesh = UnstructuredMesh(points, [Δs; □s])
```

```@example overview
# the elements can be materialized from the global vector of points
collect(elements(mesh))
```
