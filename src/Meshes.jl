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

# Queryverse compatibility
import TableTraits
import IteratorInterfaceExtensions
const IIE = IteratorInterfaceExtensions

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

# basic geometries
include("geometries.jl")

# collections
include("collections.jl")

# geometric distances
include("distances.jl")

# rotation conventions
include("conventions.jl")

# connectivities and meshes
include("connectivities.jl")
include("topostructures.jl")
include("toporelations.jl")
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
include("simplification.jl")
include("refinement.jl")
include("boundingboxes.jl")
include("supportfun.jl")

# utilities
include("utils.jl")

# plot recipes
include("plotrecipes/domain.jl")
include("plotrecipes/data.jl")
include("plotrecipes/points.jl")
include("plotrecipes/geometries.jl")
include("plotrecipes/collections.jl")
include("plotrecipes/cartesiangrids.jl")
include("plotrecipes/partitions.jl")

export 
  # points
  Point, Point1, Point2, Point3,
  Point1f, Point2f, Point3f,
  embeddim, coordtype, coordinates,
  ⪯, ≺, ⪰, ≻,

  # vectors
  Vec, Vec1, Vec2, Vec3, Vec1f, Vec2f, Vec3f,

  # angles
  ∠,

  # domain/data traits
  Domain, Data, Variable,
  domain, constructor, asarray,
  element, nelements,
  variables, name, mactype,

  # optional traits
  IsGrid, isgrid,

  # geometries
  Geometry,
  embeddim, paramdim, coordtype,
  measure, area, volume, boundary,
  centroid, isconvex,

  # primitives
  Primitive,
  Line, Ray, Plane, BezierCurve,
  Box, Ball, Sphere, Cylinder,
  ncontrols, degree, Horner, DeCasteljau,
  center, radius, height, sides,
  measure, diagonal,

  # polytopes
  Polytope, Polygon, Polyhedron,
  Segment,
  Ngon, Triangle, Quadrangle,
  Pentagon, Hexagon, Heptagon,
  Octagon, Nonagon, Decagon,
  Chain, PolyArea,
  Tetrahedron, Pyramid, Hexahedron,
  vertices, nvertices,
  windingnumber, chains, segments,
  isclosed, issimple, hasholes,
  angles, innerangles, close!, open!,
  orientation, bridge,

  # orientation algorithms
  WindingOrientation,
  TriangleOrientation,

  # multi-types
  Multi,

  # collections
  Collection,
  PointSet,
  GeometrySet,

  # rotation conventions
  RotationConvention,
  GSLIB, Leapfrog, Datamine,
  TaitBryanExtr, TaitBryanIntr,
  EulerExtr, EulerIntr,
  axesseq, orientation, angleunits,
  mainaxis, isextrinsic, rotmat,

  # connectivities
  Connectivity,
  paramdim, indices,
  connect, materialize,

  # topologies
  Topology,
  vertices, nvertices,
  faces, elements, facets,
  element, nelements,
  facet, nfacets,
  FullTopology,
  HalfEdgeTopology, HalfEdge,
  half4elem, half4vert,
  half4edge, half4pair,
  edge4pair,

  # topological relations
  TopologicalRelation,
  Boundary, Coboundary, Adjacency,

  # meshes
  Mesh,
  CartesianGrid, SimpleMesh,
  vertices, nvertices, topology,
  faces, elements, facets,
  element, nelements,
  facet, nfacets,
  spacing,

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
  BallSearch,
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
  NoIntersection,
  CrossingLines,
  OverlappingLines,
  CrossingSegments,
  MidTouchingSegments,
  CornerTouchingSegments,
  OverlappingSegments,
  OverlappingBoxes,
  FaceTouchingBoxes,
  CornerTouchingBoxes,
  intersecttype,
  hasintersect,

  # discretization
  DiscretizationMethod,
  FIST, Dehn1899,
  discretize,

  # simplification
  SimplificationMethod,
  DouglasPeucker,
  simplify,

  # refinement
  RefinementMethod,
  CatmullClark,
  refine,

  # bounding boxes
  boundingbox,

  # support function
  supportfun,

  # utililities
  signarea,
  sideof

end # module
