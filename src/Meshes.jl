# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module Meshes

using Tables
using StaticArrays
using CircularArrays
using SimpleTraits
using RecipesBase
using LinearAlgebra
using Random

using IterTools: ivec
using StatsBase: Weights
using SpecialFunctions: gamma
using Distances: Euclidean, Mahalanobis, evaluate
using ReferenceFrameRotations: angle_to_dcm
using NearestNeighbors: KDTree, BallTree, knn, inrange

import Tables
import Random
import Base: values, ==, +, -, *
import StatsBase: sample
import NearestNeighbors: MinkowskiMetric

# numerical tolerances
include("tolerances.jl")

# basic types
include("vectors.jl")
include("points.jl")
include("angles.jl")

# type traits
include("traits/domain.jl")
include("traits/data.jl")
include("traits/optional.jl")

# point sets
include("pointsets.jl")

# basic geometries
include("geometries.jl")

# geometry sets
include("geometrysets.jl")

# connectivities and meshes
include("connectivities.jl")
include("conventions.jl")
include("mesh.jl")

# domain and data views
include("views.jl")

# neighborhoods and searches
include("neighborhoods.jl")
include("neighborsearch.jl")

# partitions
include("partitions.jl")

# algorithms
include("sampling.jl")
include("partitioning.jl")
include("intersections.jl")
include("discretization.jl")
include("boundingboxes.jl")

# utilities
include("utils.jl")

# plot recipes
include("plotrecipes/domain.jl")
include("plotrecipes/data.jl")
include("plotrecipes/geometries.jl")
include("plotrecipes/pointsets.jl")
include("plotrecipes/cartesiangrids.jl")
include("plotrecipes/mesh.jl")

export 
  # points
  Point, Point1, Point2, Point3, Point1f, Point2f, Point3f,
  embeddim, coordtype, coordinates,
  ⪯, ≺, ⪰, ≻,

  # vectors
  Vec, Vec1, Vec2, Vec3, Vec1f, Vec2f, Vec3f,

  # angles
  ∠,

  # domain/data traits
  Domain, Data,
  domain, asarray,

  # optional traits
  isgrid,

  # geometries
  Geometry,
  embeddim, paramdim, coordtype,
  measure, area, volume, boundary,
  centroid,

  # primitives
  Primitive,
  Box, Ball, Sphere, Cylinder,
  center, radius, height, sides,
  measure, diagonal,

  # polytopes
  Polytope, Polygon, Polyhedron,
  Segment, Triangle, Quadrangle,
  Pyramid, Tetrahedron, Hexahedron,
  Chain, PolyArea,
  vertices, nvertices,
  faces, facets,
  windingnumber, chains,
  isclosed, issimple, hasholes,
  angles, innerangles, close!, open!,
  orientation, bridge,

  # geometry sets
  GeometrySet,

  # connectivities
  Connectivity,
  polytopetype, connect, materialize,

  # ordering conventions
  OrderingConvention,
  GMSH,
  connectivities,

  # rotation conventions
  RotationConvention,
  GSLIB, Leapfrog, Datamine,
  TaitBryanExtr, TaitBryanIntr,
  EulerExtr, EulerIntr,
  axesseq, orientation, angleunits,
  mainaxis, isextrinsic, rotmat,

  # point sets
  PointSet,
  coordinates, coordinates!,
  nelements,

  # meshes
  Mesh,
  CartesianGrid, UnstructuredMesh,
  faces, elements, spacing,
  coordinates, coordinates!,
  nelements,

  # views
  DomainView, DataView,

  # neighborhoods
  NormBall, Ellipsoid,
  metric, radius,

  # neighbordhood search
  NeighborSearchMethod,
  BoundedNeighborSearchMethod,
  NeighborhoodSearch,
  KNearestSearch,
  KBallSearch,
  BoundedSearch,
  search!, search,

  # partitions
  Partition,
  subsets, metadata,

  # sampling
  SamplingMethod,
  RegularSampling,
  UniformSampling,
  WeightedSampling,
  BallSampling,
  sample,

  # partitioning
  PartitionMethod,
  PredicatePartitionMethod,
  SPredicatePartitionMethod,
  RandomPartition,
  FractionPartition,
  BlockPartition,
  BisectPointPartition,
  BisectFractionPartition,
  BallPartition,
  PlanePartition,
  DirectionPartition,
  PredicatePartition,
  SpatialPredicatePartition,
  ProductPartition,
  HierarchicalPartition,
  partition, split, →,

  # intersections
  Intersection,
  CrossingSegments,
  MidTouchingSegments,
  CornerTouchingBoxes,
  OverlappingSegments,
  NonIntersectingSegments,
  OverlappingBoxes,
  SegmentTouchingBoxes,
  CornerTouchingBoxes,
  NonIntersectingBoxes,
  intersecttype,

  # discretization
  DiscretizationMethod,
  FIST,
  discretize,

  # bounding boxes
  boundingbox,

  # utililities
  signarea,
  sideof

end # module
