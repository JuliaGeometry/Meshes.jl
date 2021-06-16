# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module Meshes

using Tables
using StaticArrays
using SparseArrays
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

# IO utils
include("ioutils.jl")

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

# meshes
include("connectivities.jl")
include("topologies.jl")
include("toporelations.jl")
include("mesh.jl")

# mesh data
include("meshdata.jl")

# miscellaneous
include("conventions.jl")
include("paths.jl")
include("distances.jl")
include("neighborhoods.jl")
include("neighborsearch.jl")
include("supportfun.jl")
include("laplacian.jl")

# views and partitions
include("views.jl")
include("partitions.jl")

# algorithms
include("viewing.jl")
include("sampling.jl")
include("partitioning.jl")
include("intersections.jl")
include("discretization.jl")
include("simplification.jl")
include("refinement.jl")
include("smoothing.jl")
include("boundingboxes.jl")

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
  embeddim, paramdim,
  coordtype, coordinates,
  ⪯, ≺, ⪰, ≻,

  # vectors
  Vec, Vec1, Vec2, Vec3, Vec1f, Vec2f, Vec3f,

  # angles
  ∠,

  # domain traits
  Domain,
  embeddim, paramdim, coordtype,
  element, nelements,

  # data traits
  Data, Variable,
  domain, constructor, asarray,
  variables, name, mactype,

  # optional traits
  IsGrid, isgrid,

  # geometries
  Geometry,
  embeddim, paramdim, coordtype,
  measure, area, volume, boundary,
  centroid, isconvex, issimplex,

  # primitives
  Primitive,
  Line, Ray, Plane, BezierCurve,
  Box, Ball, Sphere, Cylinder,
  ncontrols, degree, Horner, DeCasteljau,
  center, radius, height, sides,
  measure, diagonal,

  # polytopes
  Polytope, Polygon, Polyhedron,
  Segment, Ngon, Triangle, Quadrangle,
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

  # point or geometry alias
  PointOrGeometry,

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

  # miscellaneous
  supportfun,
  laplacematrix,

  # connectivities
  Connectivity,
  paramdim, indices,
  connect, materialize,
  issimplex,

  # topologies
  Topology,
  vertices, nvertices,
  faces, elements, facets,
  element, nelements,
  facet, nfacets,
  FullTopology,
  GridTopology,
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
  topoconvert,
  spacing,

  # mesh data
  MeshData,
  meshdata,

  # views
  DomainView, DataView,

  # partitions
  Partition,
  subsets, metadata,

  # viewing
  viewindices,
  slice,

  # sampling
  SamplingMethod,
  DiscreteSamplingMethod,
  ContinuousSamplingMethod,
  UniformSampling,
  WeightedSampling,
  BallSampling,
  RegularSampling,
  HomogeneousSampling,
  MinDistanceSampling,
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
  QuadRefinement,
  CatmullClark,
  refine,

  # smoothing
  SmoothingMethod,
  TaubinSmoothing,
  smooth,

  # bounding boxes
  boundingbox,

  # utililities
  signarea,
  sideof

end # module
