module Meshes

using Tables
using StaticArrays
using StructArrays
using IterTools
using LinearAlgebra
using EarCut_jll

using Base: @propagate_inbounds

import Base: +, -

# basic concepts
include("vectors.jl")
include("matrices.jl")
include("points.jl")

include("fixed_arrays.jl")
include("offsetintegers.jl")
include("basictypes.jl")

include("primitives/rectangles.jl")
include("primitives/spheres.jl")
include("primitives/cylinders.jl")
include("primitives/pyramids.jl")

include("interfaces.jl")
include("metadata.jl")
include("viewtypes.jl")
include("primitives.jl")
include("meshes.jl")
include("triangulation.jl")
include("lines.jl")
include("boundingboxes.jl")

# points
export Point, Point2, Point3, Point2f, Point3f

# TODO: review these
export PointMeta

# vectors
export Vec, Vec2, Vec3, Vec2f, Vec3f
export vunit, vfill

# geometries
export AbstractGeometry, GeometryPrimitive
export LineFace, Polytope, Line, NgonFace
export LineString, AbstractPolygon, Polygon, MultiPoint, MultiLineString, MultiPolygon
export Simplex, connect, Triangle, NSimplex, Tetrahedron
export QuadFace, metafree, coordinates, TetrahedronFace
export TupleView, SimplexFace, Mesh, meta
export Triangle, TriangleP
export AbstractFace, TriangleFace, QuadFace, GLTriangleFace
export OffsetInteger, ZeroIndex, OneIndex, GLIndex
export FaceView, SimpleFaceView
export PolygonMeta, MultiPointMeta, MultiLineStringMeta, MeshMeta, LineStringMeta,
       MultiPolygonMeta
export decompose, coordinates, faces, normals, decompose_uv, decompose_normals,
       texturecoordinates, convert_simplex
export Tesselation, pointmeta, Normal, UV, UVW
export GLTriangleFace, GLNormalMesh3D, GLPlainTriangleMesh, GLUVMesh3D, GLUVNormalMesh3D
export AbstractMesh, Mesh, TriangleMesh
export GLNormalMesh2D, PlainTriangleMesh
export MetaT, meta_table

# primitives
export Rectangle, Cylinder, Pyramid, Sphere

# TODO: review these
export HyperSphere, Circle
export Cylinder2, Cylinder3

export height, origin, radius, width, widths, xwidth, yheight
export before, during, isinside, isoutside, meets, overlaps, intersects, finishes
export centered, direction, area, update, extremity
export max_dist_dim, max_euclidean, max_euclideansq, min_dist_dim, min_euclidean
export min_euclideansq, minmax_dist_dim, minmax_euclidean, minmax_euclideansq
export self_intersections, split_intersections

# meshes
export AbstractMesh, TriangleMesh, PlainMesh, GLPlainMesh, GLPlainMesh2D, GLPlainMesh3D
export UVMesh, GLUVMesh, GLUVMesh2D, GLUVMesh3D
export NormalMesh, GLNormalMesh, GLNormalMesh2D, GLNormalMesh3D
export NormalUVMesh, GLNormalUVMesh, GLNormalUVMesh2D, GLNormalUVMesh3D
export NormalUVWMesh, GLNormalUVWMesh, GLNormalUVWMesh2D, GLNormalUVWMesh3D

# mesh creation
export triangle_mesh, triangle_mesh, uv_mesh
export uv_mesh, normal_mesh, uv_normal_mesh

# bounding boxes
export boundingbox

end # module
