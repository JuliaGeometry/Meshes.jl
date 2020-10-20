# Meshes.jl

[![Build Status](https://travis-ci.com/JuliaGeometry/Meshes.jl.svg?branch=master)](https://travis-ci.com/JuliaGeometry/Meshes.jl)
[![Codecov](https://codecov.io/gh/JuliaGeometry/Meshes.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGeometry/Meshes.jl)

This package provides efficient implementations of concepts from
mathematical geometry and finite element analysis, which are useful
in geospatial modeling and numerical simulation. It was forked from
the GeometryBasics.jl package after careful considerations.

Unlike GeometryBasics.jl, which was originally designed with visualization
pipelines in mind, Meshes.jl is concerned with rigorous mathematical
descriptions of geometries and finite element discretizations (a.k.a.
meshes) that are adequate for general manifolds embedded in Rⁿ including
surfaces described with spherical coordinates, and geometries described
with different coordinate reference systems.

For more background on the necessity of this fork, please check some of
the following discussions:

- Point vs. Vector: https://github.com/JuliaGeometry/GeometryBasics.jl/issues/91
- Polygon as vector: https://github.com/JuliaGeometry/GeometryBasics.jl/issues/95
- Mesh as vector: https://github.com/JuliaGeometry/GeometryBasics.jl/issues/54
- Lack of supertype: https://github.com/JuliaGeometry/GeometryBasics.jl/issues/70
- Lack of interface: https://github.com/JuliaGeometry/GeometryBasics.jl/issues/15

## Usage

The code is being completely refactored to accommodate advanced use cases.
Until then, we kindly ask users to read the source code, the test suite,
and the docstrings that we are adding as part of this effort.

## Contributing

Contributions are very welcome, as are feature requests and suggestions.
Please open an issue if you encounter any problems. We take issues very
seriously and value any type of feedback.

We have high standards when it comes to source code. Please adopt the
coding style that is present in the files when submitting pull requests.
