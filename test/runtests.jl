using Meshes
using Test, Pkg, Random

# workaround GR warnings
ENV["GKSwstype"] = "100"

# environment settings
islinux = Sys.islinux()
istravis = "TRAVIS" ∈ keys(ENV)
isappveyor = "APPVEYOR" ∈ keys(ENV)
isCI = istravis || isappveyor
datadir = joinpath(@__DIR__,"data")

# list of tests
testfiles = [
  "points.jl",
  "angles.jl",
  "geometries.jl",
  "polytopes.jl",
  "primitives.jl",
  "mesh.jl",
  "sampling.jl",
  "discretization.jl",
  "boundingboxes.jl"
]

# --------------------------------
# RUN TESTS WITH SINGLE PRECISION
# --------------------------------
T = Float32
P1, P2, P3 = Point{2,T}, Point{3,T}
V2, V3 = Vec{2,T}, Vec{3,T}
@testset "Meshes.jl ($T)" begin
  for testfile in testfiles
    include(testfile)
  end
end

# --------------------------------
# RUN TESTS WITH DOUBLE PRECISION
# --------------------------------
T = Float64
P1, P2, P3 = Point{2,T}, Point{3,T}
V2, V3 = Vec{2,T}, Vec{3,T}
@testset "Meshes.jl ($T)" begin
  for testfile in testfiles
    include(testfile)
  end
end
