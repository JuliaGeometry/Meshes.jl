# Meshes.jl

[![Build Status](https://travis-ci.com/JuliaGeometry/Meshes.jl.svg?branch=master)](https://travis-ci.com/JuliaGeometry/Meshes.jl)
[![Codecov](https://codecov.io/gh/JuliaGeometry/Meshes.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaGeometry/Meshes.jl)

This project provides efficient implementations of concepts from
computational geometry and finite element analysis. It promotes
rigorous mathematical definitions of spatial discretizations
(a.k.a. meshes) that are adequate for describing general
manifolds embedded in Rⁿ, including surfaces described with
spherical coordinates, and geometries described with multiple
coordinate reference systems.

Our long-term goal is to provide all the features of the
[CGAL](https://www.cgal.org) project in pure Julia.

## Usage

The project is being actively developed. We kindly ask users to
read the source code, the tests, and the docstrings for usage
examples.

## Contributing

Contributions are very welcome, as are feature requests and suggestions.
Please open an issue if you encounter any problems. We take issues very
seriously and value any type of feedback.

We have high standards when it comes to source code. Please adopt the
coding style that is present in the files when submitting pull requests.

### Related packages

This project is an evolution of previous efforts:

- [GeometryTypes.jl](https://github.com/JuliaGeometry/GeometryTypes.jl)
- [GeometryBasics.jl](https://github.com/JuliaGeometry/GeometryBasics.jl)

The main difference is that we provide robust algorithms that work well
in the context of finite element analysis and geospatial modeling. 
