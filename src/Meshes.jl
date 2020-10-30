# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module Meshes

using IterTools
using StaticArrays
using LinearAlgebra

import Base: ==, +, -

# basic types
include("vectors.jl")
include("matrices.jl")
include("points.jl")

# geometries and meshes
include("geometries.jl")
include("connectivities.jl")
include("meshes.jl")

# algorithms
include("boundingboxes.jl")

export 
  # points
  Point, Point2, Point3, Point2f, Point3f,
  embeddim, coordtype, coordinates,

  # vectors
  Vec, Vec2, Vec3, Vec2f, Vec3f,

  # geometries
  Geometry,
  embeddim, coordtype, boundary,

  # primitives
  Primitive,
  Box, Ball, Sphere, Cylinder,
  center, radius, height, sides, volume,

  # polytopes
  Polytope,
  vertices,
  facets,

  # faces
  Face,
  Segment, Triangle, Quadrangle,
  Pyramid, Tetrahedron, Hexahedron,
  rank,

  # chains
  Chain,
  isclosed,

  # polygons
  Polygon,

  # sets
  GeometrySet,

  # connectivities
  Connectivity,
  facetype, connect, materialize,

  # meshes
  Mesh,
  CartesianGrid, UnstructuredMesh,
  faces, elements, spacing,

  # bounding boxes
  boundingbox

end # module
