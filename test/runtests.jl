using Meshes
using Tables
using Distances
using Statistics
using LinearAlgebra
using CoordRefSystems
using CategoricalArrays
using CircularArrays
using StaticArrays
using SparseArrays
using PlyIO
using Unitful
using Rotations
using Test, StableRNGs
using ReferenceTests, ImageIO

using TransformsBase: Identity, →

import TransformsBase as TB
import CairoMakie as Mke

# environment settings
isCI = "CI" ∈ keys(ENV)
islinux = Sys.islinux()
visualtests = !isCI || (isCI && islinux)
datadir = joinpath(@__DIR__, "data")

# dummy definitions
include("dummy.jl")

# helper functions
include("testutils.jl")

cart(args...) = cart(T, args...)

latlon(args...) = latlon(T, args...)

vector(args...) = vector(T, args...)

cartgrid(args...) = cartgrid(T, args...)

randpoint1(n) = randcart(T, 1, n)
randpoint2(n) = randcart(T, 2, n)
randpoint3(n) = randcart(T, 3, n)

# list of tests
testfiles = [
  "vectors.jl",
  "primitives.jl",
  "polytopes.jl",
  "multigeoms.jl",
  "connectivities.jl",
  "topologies.jl",
  "toporelations.jl",
  "domains.jl",
  "subdomains.jl",
  "sets.jl",
  "mesh.jl",
  "trajecs.jl",
  "utils.jl",
  "viewing.jl",
  "partitioning.jl",
  "sorting.jl",
  "traversing.jl",
  "neighborhoods.jl",
  "neighborsearch.jl",
  "predicates.jl",
  "winding.jl",
  "sideof.jl",
  "orientation.jl",
  "merging.jl",
  "clipping.jl",
  "clamping.jl",
  "intersections.jl",
  "complement.jl",
  "simplification.jl",
  "boundingboxes.jl",
  "hulls.jl",
  "sampling.jl",
  "pointification.jl",
  "tesselation.jl",
  "discretization.jl",
  "refinement.jl",
  "coarsening.jl",
  "transforms.jl",
  "distances.jl",
  "supportfun.jl",
  "matrices.jl",
  "tolerances.jl"
]

# --------------------------------
# RUN TESTS WITH SINGLE PRECISION
# --------------------------------
T = Float32
ℳ = Meshes.Met{T}
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
ℳ = Meshes.Met{T}
@testset "Meshes.jl ($T)" begin
  for testfile in testfiles
    println("Testing $testfile...")
    include(testfile)
  end
end
