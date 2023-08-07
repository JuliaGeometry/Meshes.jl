# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module Meshes

using Tables
using StaticArrays
using SparseArrays
using CircularArrays
using LinearAlgebra
using Random

using Bessels: gamma
using StatsBase: AbstractWeights, Weights, quantile
using Distances: PreMetric, Euclidean, Mahalanobis, evaluate
using Rotations: Rotation, QuatRotation, Angle2d, rotation_between
using NearestNeighbors: KDTree, BallTree, knn, inrange

import Tables
import Random
import Base: values
import Base: ==, !
import Base: +, -, *
import StatsBase: sample
import Distances: evaluate
import NearestNeighbors: MinkowskiMetric

# Transforms API
import TransformsBase: Transform, →
import TransformsBase: isrevertible, isinvertible
import TransformsBase: apply, revert, reapply

# IO utils
include("ioutils.jl")

# numerical tolerances
include("tolerances.jl")

# basic vector type
include("vectors.jl")

# geometries
include("geometries.jl")

# topological objects
include("connectivities.jl")
include("topologies.jl")
include("toporelations.jl")

# domains
include("domains.jl")

# data over domains
include("data.jl")

# utilities
include("utils.jl")

# domain views
include("views.jl")
include("viewing.jl")

# domain partitions
include("partitions.jl")
include("partitioning.jl")

# domain traversal
include("traversing.jl")

# neighbor search
include("neighborhoods.jl")
include("neighborsearch.jl")

# miscellaneous
include("distances.jl")
include("supportfun.jl")
include("matrices.jl")
include("projecting.jl")

# predicates
include("predicates.jl")

# operations
include("merging.jl")
include("clipping.jl")
include("intersections.jl")
include("complement.jl")
include("simplification.jl")
include("boundingboxes.jl")
include("hulls.jl")
include("sampling.jl")
include("pointification.jl")
include("discretization.jl")
include("refinement.jl")

# transforms
include("transforms.jl")

# visualization
include("viz.jl")

export
  # vectors
  Vec,
  Vec1,
  Vec2,
  Vec3,
  Vec1f,
  Vec2f,
  Vec3f,
  ∠,
  ⋅,
  ×,

  # geometries
  Geometry,
  embeddim,
  paramdim,
  coordtype,
  measure,
  area,
  volume,
  boundary,
  center,
  centroid,
  perimeter,

  # primitives
  Primitive,
  Point,
  Point1,
  Point2,
  Point3,
  Point1f,
  Point2f,
  Point3f,
  Ray,
  Line,
  BezierCurve,
  Plane,
  Box,
  Ball,
  Sphere,
  Disk,
  Circle,
  Cylinder,
  CylinderSurface,
  Cone,
  ConeSurface,
  Torus,
  controls,
  ncontrols,
  degree,
  Horner,
  DeCasteljau,
  coordinates,
  radius,
  radii,
  plane,
  bottom,
  top,
  axis,
  isright,
  sides,
  measure,
  diagonal,
  ⪯,
  ≺,
  ⪰,
  ≻,

  # polytopes
  Polytope,
  Chain,
  Segment,
  Rope,
  Ring,
  Polygon,
  Ngon,
  Triangle,
  Quadrangle,
  Pentagon,
  Hexagon,
  Heptagon,
  Octagon,
  Nonagon,
  Decagon,
  PolyArea,
  Polyhedron,
  Tetrahedron,
  Pyramid,
  Hexahedron,
  vertex,
  vertices,
  nvertices,
  windingnumber,
  rings,
  segments,
  angles,
  innerangles,
  orientation,
  bridge,
  normal,

  # multi-geometries
  Multi,
  MultiPoint,
  MultiSegment,
  MultiRope,
  MultiRing,
  MultiPolygon,

  # connectivities
  Connectivity,
  paramdim,
  indices,
  connect,
  materialize,
  pltype,

  # topologies
  Topology,
  vertex,
  vertices,
  nvertices,
  element,
  elements,
  nelements,
  facet,
  facets,
  nfacets,
  faces,
  nfaces,
  GridTopology,
  HalfEdgeTopology,
  HalfEdge,
  SimpleTopology,
  elem2cart,
  cart2elem,
  corner2cart,
  cart2corner,
  elem2corner,
  corner2elem,
  elementtype,
  facettype,
  half4elem,
  half4vert,
  half4edge,
  half4pair,
  edge4pair,
  connec4elem,

  # topological relations
  TopologicalRelation,
  Boundary,
  Coboundary,
  Adjacency,

  # domain traits
  Domain,
  embeddim,
  paramdim,
  coordtype,
  element,
  nelements,

  # sets
  GeometrySet,
  PointSet,

  # meshes
  Mesh,
  Grid,
  CartesianGrid,
  RectilinearGrid,
  SimpleMesh,
  vertex,
  vertices,
  nvertices,
  element,
  elements,
  nelements,
  facet,
  facets,
  nfacets,
  topology,
  topoconvert,
  cart2vert,
  spacing,
  offset,

  # data traits
  Data,
  domain,
  constructor,
  asarray,

  # mesh data
  MeshData,
  meshdata,

  # domain/data traits
  nitems,

  # views
  DomainView,
  DataView,

  # viewing
  unview,
  indices,

  # partitions
  Partition,
  indices,
  metadata,

  # partitioning
  PartitionMethod,
  PredicatePartitionMethod,
  SPredicatePartitionMethod,
  UniformPartition,
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
  partition,
  split,

  # traversing
  Path,
  LinearPath,
  RandomPath,
  ShiftedPath,
  SourcePath,
  MultiGridPath,
  traverse,

  # neighborhoods
  Neighborhood,
  MetricBall,
  metric,
  radii,
  radius,
  isisotropic,

  # neighborhood search
  NeighborSearchMethod,
  BoundedNeighborSearchMethod,
  BallSearch,
  KNearestSearch,
  KBallSearch,
  BoundedSearch,
  GlobalSearch,
  search!,
  search,
  maxneighbors,

  # predicates
  isparametrized,
  isperiodic,
  issimplex,
  isclosed,
  isconvex,
  issimple,
  hasholes,
  intersects,
  iscollinear,
  iscoplanar,

  # clipping
  ClippingMethod,
  SutherlandHodgman,
  clip,

  # intersections
  IntersectionType,
  Crossing,
  CornerCrossing,
  EdgeCrossing,
  Touching,
  CornerTouching,
  EdgeTouching,
  Overlapping,
  PosOverlapping,
  NegOverlapping,
  NotIntersecting,
  Intersecting,

  # intersecting
  Intersection,
  intersection,
  type,

  # simplification
  SimplificationMethod,
  DouglasPeucker,
  Selinger,
  simplify,
  decimate,

  # bounding boxes
  boundingbox,

  # hulls
  HullMethod,
  GrahamScan,
  JarvisMarch,
  hull,
  convexhull,

  # sampling
  SamplingMethod,
  DiscreteSamplingMethod,
  ContinuousSamplingMethod,
  UniformSampling,
  WeightedSampling,
  BallSampling,
  BlockSampling,
  RegularSampling,
  HomogeneousSampling,
  MinDistanceSampling,
  sampleinds,
  sample,

  # pointification
  pointify,

  # discretization
  DiscretizationMethod,
  BoundaryDiscretizationMethod,
  FanTriangulation,
  RegularDiscretization,
  FIST,
  Dehn1899,
  Tetrahedralization,
  discretize,
  discretizewithin,
  simplexify,

  # refinement
  RefinementMethod,
  TriRefinement,
  QuadRefinement,
  CatmullClark,
  TriSubdivision,
  refine,

  # transforms
  GeometricTransform,
  CoordinateTransform,
  Rotate,
  Translate,
  Stretch,
  StdCoords,
  Repair,
  LambdaMuSmoothing,
  LaplaceSmoothing,
  TaubinSmoothing,

  # miscellaneous
  WindingOrientation,
  TriangleOrientation,
  signarea,
  sideof,
  householderbasis,
  supportfun,
  laplacematrix,
  measurematrix,
  adjacencymatrix,
  atol,

  # visualization
  viz,
  viz!

end # module
