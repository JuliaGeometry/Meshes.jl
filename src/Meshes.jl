# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module Meshes

using StaticArrays
using IterTools
using LinearAlgebra
using EarCut_jll

using Base: @propagate_inbounds

import Base: +, -

# basic concepts
include("vectors.jl")
include("matrices.jl")
include("points.jl")

include("geometries.jl")

include("interfaces.jl")
include("viewtypes.jl")
include("primitives.jl")
include("meshes.jl")
include("triangulation.jl")
include("lines.jl")
include("boundingboxes.jl")

# points
export Point, Point2, Point3, Point2f, Point3f
export coordtype, coordinates

# vectors
export Vec, Vec2, Vec3, Vec2f, Vec3f
export vunit, vfill

# geometries
export Geometry, Primitive, Polytope
export coordtype

# primitives
export Box, Sphere, Cylinder
export center, radius, height, sides, volume

export LineFace, Line, NgonFace
export LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon
export Simplex, connect, Triangle, NSimplex, Tetrahedron
export AbstractFace, TriangleFace, QuadFace, TetrahedronFace
export FaceView, SimpleFaceView, TupleView
export decompose, faces, convert_simplex
export Tesselation

export before, during, isinside, isoutside, overlaps, intersects, finishes
export direction, area, update, extremity
export self_intersections, split_intersections

# meshes
export AbstractMesh, Mesh, TriangleMesh, PlainMesh, PlainTriangleMesh

# mesh creation
export triangle_mesh

# bounding boxes
export boundingbox

end # module
