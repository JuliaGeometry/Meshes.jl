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
    vunit, vfill,

    # geometries
    Geometry,
    coordtype,

    # primitives
    Primitive,
    Box, Ball, Sphere, Cylinder,
    center, radius, height, sides, volume,

    # polytopes
    Polytope,
    Segment, Triangle, Quadrangle,
    Pyramid, Tetrahedron, Hexahedron,
    vertices,

    # chains
    Chain,

    # connectivities
    Connectivity,

    # meshes
    Mesh,
    UnstructuredMesh,

    # TODO: review these
    before, during, overlaps, finishes,

    # bounding boxes
    boundingbox

end # module
