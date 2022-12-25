using Meshes
using Tables
using Distances
using Statistics
using LinearAlgebra
using CategoricalArrays
using CircularArrays
using StaticArrays
using SparseArrays
using PlyIO
using Unitful
using MeshViz
using ReferenceFrameRotations
using Test, Random
using ReferenceTests, ImageIO

import TransformsBase as TB
import CairoMakie as Mke

# environment settings
isCI = "CI" ∈ keys(ENV)
islinux = Sys.islinux()
visualtests = !isCI || (isCI && islinux)
datadir = joinpath(@__DIR__,"data")

# helper function to read *.line files containing polygons
# generated with RPG (https://github.com/cgalab/genpoly-rpg)
function readpoly(T, fname)
  open(fname, "r") do f
    # read outer chain
    n = parse(Int, readline(f))
    outer = map(1:n) do _
      coords = readline(f)
      x, y = parse.(T, split(coords))
      Point(x, y)
    end

    # read inner chains
    inners = []
    while !eof(f)
      n = parse(Int, readline(f))
      inner = map(1:n) do _
        coords = readline(f)
        x, y = parse.(T, split(coords))
        Point(x, y)
      end
      push!(inners, inner)
    end

    # return polygonal area
    PolyArea(outer, inners)
  end
end

# helper function to read *.ply files containing meshes
function readply(T, fname)
  ply = load_ply(fname)
  x = ply["vertex"]["x"]
  y = ply["vertex"]["y"]
  z = ply["vertex"]["z"]
  points = Point{3,T}.(x, y, z)
  connec = [connect(Tuple(c.+1)) for c in ply["face"]["vertex_indices"]]
  SimpleMesh(points, connec)
end

include("dummy.jl")

# list of tests
testfiles = [
  "vectors.jl",
  "points.jl",
  "angles.jl",
  "polytopes.jl",
  "primitives.jl",
  "multi.jl",
  "collections.jl",
  "connectivities.jl",
  "topologies.jl",
  "toporelations.jl",
  "mesh.jl",
  "meshdata.jl",
  "traits.jl",
  "paths.jl",
  "distances.jl",
  "neighborhoods.jl",
  "neighborsearch.jl",
  "supportfun.jl",
  "matrices.jl",
  "views.jl",
  "viewing.jl",
  "merging.jl",
  "sampling.jl",
  "partitioning.jl",
  "intersections.jl",
  "discretization.jl",
  "simplification.jl",
  "refinement.jl",
  "boundingboxes.jl",
  "hulls.jl",
  "transforms.jl",
  "utils.jl",
  "cleaning.jl"
]

# --------------------------------
# RUN TESTS WITH SINGLE PRECISION
# --------------------------------
T = Float32
P1, P2, P3 = Point{1,T}, Point{2,T}, Point{3,T}
V1, V2, V3 = Vec{1,T}, Vec{2,T}, Vec{3,T}
@testset "Meshes.jl ($T)" begin
  for testfile in testfiles
    println("Testing $testfile...")
    include(testfile)
  end
end

# --------------------------------
# RUN TESTS WITH DOUBLE PRECISION
# --------------------------------
T = Float64
P1, P2, P3 = Point{1,T}, Point{2,T}, Point{3,T}
V1, V2, V3 = Vec{1,T}, Vec{2,T}, Vec{3,T}
@testset "Meshes.jl ($T)" begin
  for testfile in testfiles
    println("Testing $testfile...")
    include(testfile)
  end
end
