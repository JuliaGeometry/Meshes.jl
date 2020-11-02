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

# basic types
include("vectors.jl")
include("points.jl")

# geometries and meshes
include("geometries.jl")
include("connectivities.jl")
include("meshes.jl")

# algorithms
include("sampling.jl")
include("boundingboxes.jl")

# plot recipes
include("plotrecipes/chains.jl")
include("plotrecipes/polygons.jl")

export 
  # points
  Point, Point2, Point3, Point2f, Point3f,
  embeddim, coordtype, coordinates,

  # vectors
  Vec, Vec2, Vec3, Vec2f, Vec3f,

  # geometries
  Geometry,
  embeddim, paramdim, coordtype,
  measure, boundary,

  # primitives
  Primitive,
  Box, Ball, Sphere, Cylinder,
  center, radius, height, sides,

  # polytopes
  Polytope,
  vertices, nvertices,
  facets,

  # faces
  Face,
  Segment, Triangle, Quadrangle,
  Pyramid, Tetrahedron, Hexahedron,

  # chains
  Chain,
  isclosed,

  # polygons
  Polygon,
  rings, hasholes,

  # sets
  GeometrySet,

  # connectivities
  Connectivity,
  facetype, connect, materialize,

  # meshes
  Mesh,
  CartesianGrid, UnstructuredMesh,
  faces, elements, spacing,

  # sampling
  SamplingMethod,
  RegularSampling,
  sample,

  # bounding boxes
  boundingbox

end # module
