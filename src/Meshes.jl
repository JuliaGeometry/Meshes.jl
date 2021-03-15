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
using Distances: PreMetric, Euclidean, Mahalanobis, evaluate
using ReferenceFrameRotations: angle_to_dcm
using NearestNeighbors: KDTree, BallTree, knn, inrange

# import categorical arrays as a temporary solution for plot recipes
using CategoricalArrays: CategoricalValue, levelcode

import Tables
import Random
import Base: values, ==, +, -, *
import StatsBase: sample
import Distances: evaluate
import NearestNeighbors: MinkowskiMetric

# numerical tolerances
include("tolerances.jl")

# basic types
include("vectors.jl")
include("points.jl")
include("angles.jl")

# type traits
include("traits/variable.jl")
include("traits/domain.jl")
include("traits/data.jl")
include("traits/optional.jl")

# point sets
include("pointsets.jl")

# basic geometries
include("geometries.jl")

# geometric distances
include("distances.jl")

# geometry sets
include("geometrysets.jl")

# connectivities and meshes
include("connectivities.jl")
include("conventions.jl")
include("mesh.jl")

# domain and data views
include("views.jl")

# paths for traversal
include("paths.jl")

# neighborhoods and searches
include("neighborhoods.jl")
include("neighborsearch.jl")

# partitions
include("partitions.jl")

# algorithms
include("viewing.jl")
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
include("plotrecipes/partitions.jl")

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
  Domain, Data, Variable,
  domain, constructor, asarray,
  variables, name, mactype,

  # optional traits
  IsGrid, isgrid,

  # geometries
  Geometry,
  embeddim, paramdim, coordtype,
  measure, area, volume, boundary,
  centroid,

  # primitives
  Primitive,
  Line, Box, Ball, Sphere, Cylinder,
  center, radius, height, sides,
  points, measure, diagonal,

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

  # orientation algorithms
  WindingOrientation,
  TriangleOrientation,

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
  nelements,

  # meshes
  Mesh,
  CartesianGrid, UnstructuredMesh,
  faces, elements, spacing,
  nelements,

  # views
  DomainView, DataView,

  # paths
  Path,
  LinearPath, RandomPath,
  ShiftedPath, SourcePath,
  traverse,

  # neighborhoods
  Neighborhood, MetricBall,
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

  # viewing
  viewindices,
  slice,

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
  FaceTouchingBoxes,
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
