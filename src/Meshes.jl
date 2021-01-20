# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module Meshes

using IterTools
using StaticArrays
using SpecialFunctions
using LinearAlgebra
using RecipesBase

import Base: ==, +, -
import Random

# numerical tolerances
include("tolerances.jl")

# basic types
include("vectors.jl")
include("points.jl")
include("angles.jl")

# geometries and meshes
include("geometries.jl")
include("connectivities.jl")
include("conventions.jl")
include("mesh.jl")

# algorithms
include("sampling.jl")
include("discretization.jl")
include("boundingboxes.jl")

# plot recipes
include("plotrecipes/polytopes.jl")
include("plotrecipes/chains.jl")
include("plotrecipes/polyareas.jl")

export 
  # points
  Point, Point2, Point3, Point2f, Point3f,
  embeddim, coordtype, coordinates,

  # vectors
  Vec, Vec2, Vec3, Vec2f, Vec3f,

  # angles
  ∠,

  # geometries
  Geometry,
  GeometrySet,
  embeddim, paramdim, coordtype,
  measure, boundary,

  # primitives
  Primitive,
  Box, Ball, Sphere, Cylinder,
  center, radius, height, sides,

  # polytopes
  Polytope, Polygon, Polyhedron,
  Segment, Triangle, Quadrangle,
  Pyramid, Tetrahedron, Hexahedron,
  Chain, PolyArea,
  vertices, nvertices,
  faces, facets, center,
  windingnumber, chains,
  isclosed, issimple, hasholes,
  angles, innerangles, close!, open!,
  orientation, bridge,

  # connectivities
  Connectivity,
  polytopetype, connect, materialize,

  # conventions
  OrderingConvention,
  GMSH,
  connectivities,

  # meshes
  Mesh,
  CartesianGrid, UnstructuredMesh,
  faces, elements, spacing,

  # sampling
  SamplingMethod,
  RegularSampling,
  sample,

  # discretization
  DiscretizationMethod,
  FIST,
  discretize,

  # bounding boxes
  boundingbox

end # module
