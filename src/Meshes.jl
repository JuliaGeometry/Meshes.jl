# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module Meshes

using StaticArrays
using IterTools
using LinearAlgebra

import Base: +, -

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
  coordtype, coordinates,

  # vectors
  Vec, Vec2, Vec3, Vec2f, Vec3f,

  # geometries
  Geometry,
  coordtype,

  # primitives
  Primitive,
  Box, Ball, Sphere, Cylinder,
  center, radius, height, sides, volume,

  # chains
  Chain,
  vertices, isclosed,

  # polytopes
  Polytope,
  vertices,

  # faces
  Face,
  Segment, Triangle, Quadrangle,
  Pyramid, Tetrahedron, Hexahedron,
  vertices, rank,

  # polygons
  Polygon,

  # connectivities
  Connectivity,
  connect, materialize,

  # meshes
  Mesh,
  UnstructuredMesh,

  # TODO: review these
  before, during, overlaps, finishes,

  # bounding boxes
  boundingbox

end # module
