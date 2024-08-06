# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module Meshes

using CoordRefSystems
using StaticArrays
using SparseArrays
using CircularArrays
using LinearAlgebra
using Unitful
using Random

using Bessels: gamma
using Unitful: AbstractQuantity, numtype
using StatsBase: AbstractWeights, Weights, quantile
using Distances: PreMetric, Euclidean, Mahalanobis
using Distances: Haversine, SphericalAngle
using Distances: evaluate, result_type
using Rotations: Rotation, QuatRotation, Angle2d
using Rotations: rotation_between
using NearestNeighbors: KDTree, BallTree
using NearestNeighbors: knn, inrange
using DelaunayTriangulation: triangulate, voronoi
using DelaunayTriangulation: each_solid_triangle
using DelaunayTriangulation: get_polygons
using DelaunayTriangulation: get_polygon_points
using ScopedValues: ScopedValue
using Base.Cartesian: @nloops, @nref, @ntuple
using Base: @propagate_inbounds

import Random
import Base: sort
import Base: ==, !
import Base: +, -, *
import Base: <, >, ≤, ≥
import StatsBase: sample
import Distances: evaluate
import NearestNeighbors: MinkowskiMetric

# Transforms API
import TransformsBase: Transform, →
import TransformsBase: isrevertible, isinvertible
import TransformsBase: apply, revert, reapply, inverse
import TransformsBase: parameters, preprocess

# CoordRefSystems API
import CoordRefSystems: lentype

# unit utils
include("units.jl")

# IO utils
include("ioutils.jl")

# numerical tolerances
include("tolerances.jl")

# basic vector type
include("vectors.jl")

# manifold types
include("manifolds.jl")

# geometries
include("geometries.jl")

# topological objects
include("connectivities.jl")
include("topologies.jl")
include("toporelations.jl")

# domains
include("domains.jl")

# utilities
include("utils.jl")

# domain views
include("viewing.jl")

# domain partitions
include("partitions.jl")
include("partitioning.jl")

# domain sorting
include("sorting.jl")

# domain traversal
include("traversing.jl")

# neighbor search
include("neighborhoods.jl")
include("neighborsearch.jl")

# predicates
include("predicates.jl")

# operations
include("winding.jl")
include("sideof.jl")
include("orientation.jl")
include("measures.jl")
include("boundary.jl")
include("merging.jl")
include("clipping.jl")
include("clamping.jl")
include("intersections.jl")
include("complement.jl")
include("simplification.jl")
include("boundingboxes.jl")
include("hulls.jl")
include("sampling.jl")
include("pointification.jl")
include("tesselation.jl")
include("discretization.jl")
include("refinement.jl")
include("coarsening.jl")

# transforms
include("transforms.jl")

# miscellaneous
include("rand.jl")
include("distances.jl")
include("supportfun.jl")
include("matrices.jl")
include("projecting.jl")

# visualization
include("viz.jl")

export
  # vectors
  Vec,
  ∠,
  ⋅,
  ×,

  # manifolds
  𝔼,
  🌐,

  # geometries
  Geometry,
  embeddim,
  paramdim,
  crs,
  manifold,
  center,
  centroid,

  # primitives
  Primitive,
  Point,
  Ray,
  Line,
  BezierCurve,
  Plane,
  Box,
  Ball,
  Sphere,
  Ellipsoid,
  Disk,
  Circle,
  Cylinder,
  CylinderSurface,
  ParaboloidSurface,
  Cone,
  ConeSurface,
  Frustum,
  FrustumSurface,
  Torus,
  controls,
  ncontrols,
  degree,
  Horner,
  DeCasteljau,
  coords,
  to,
  radius,
  radii,
  plane,
  bottom,
  top,
  axis,
  base,
  apex,
  isright,
  sides,
  diagonal,
  focallength,

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
  Hexahedron,
  Pyramid,
  Wedge,
  vertex,
  vertices,
  nvertices,
  rings,
  segments,
  angles,
  innerangles,
  normal,
  ≗,

  # multi-geometries
  Multi,
  MultiPoint,
  MultiSegment,
  MultiRope,
  MultiRing,
  MultiPolygon,
  MultiPolyhedron,

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
  SubDomain,
  embeddim,
  paramdim,
  crs,
  manifold,
  element,
  nelements,

  # sets
  GeometrySet,
  PointSet,

  # meshes
  Mesh,
  Grid,
  SubGrid,
  CartesianGrid,
  RectilinearGrid,
  StructuredGrid,
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
  spacing,
  offset,

  # trajectories
  CylindricalTrajectory,

  # viewing
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

  # sorting
  SortingMethod,
  DirectionSort,
  sortinds,

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
  radii,
  rotation,
  metric,
  radius,
  isisotropic,

  # neighborhood search
  NeighborSearchMethod,
  BoundedNeighborSearchMethod,
  BallSearch,
  KNearestSearch,
  KBallSearch,
  search!,
  searchdists!,
  search,
  searchdists,
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
  ≺,
  ≻,
  ⪯,
  ⪰,

  # winding number
  winding,

  # sideof
  sideof,
  SideType,
  IN,
  OUT,
  ON,
  LEFT,
  RIGHT,

  # orientation
  orientation,
  OrientationType,
  CW,
  CCW,

  # measures
  measure,
  area,
  volume,
  perimeter,

  # boundary
  boundary,

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

  # tesselation
  TesselationMethod,
  DelaunayTesselation,
  VoronoiTesselation,
  tesselate,

  # discretization
  DiscretizationMethod,
  BoundaryTriangulationMethod,
  FanTriangulation,
  DehnTriangulation,
  HeldTriangulation,
  DelaunayTriangulation,
  ManualDiscretization,
  RegularDiscretization,
  discretize,
  discretizewithin,
  simplexify,

  # refinement
  RefinementMethod,
  TriRefinement,
  QuadRefinement,
  RegularRefinement,
  CatmullClark,
  TriSubdivision,
  refine,

  # coarsening
  CoarseningMethod,
  RegularCoarsening,
  coarsen,

  # transforms
  GeometricTransform,
  CoordinateTransform,
  Rotate,
  Translate,
  Scale,
  Affine,
  Stretch,
  StdCoords,
  Proj,
  LengthUnit,
  Shadow,
  Within,
  Repair,
  Bridge,
  LambdaMuSmoothing,
  LaplaceSmoothing,
  TaubinSmoothing,
  isaffine,
  isinvertible,
  inverse,

  # miscellaneous
  signarea,
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
