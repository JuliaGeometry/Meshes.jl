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
  "faces.jl",
  "polygons.jl",
  "primitives.jl",
  "meshes.jl",
  "sampling.jl",
  "boundingboxes.jl"
]

@testset "Meshes.jl" begin
  for testfile in testfiles
    # run with single precision
    global T = Float32
    global P2 = Point{2,T}
    global P3 = Point{3,T}
    include(testfile)

    # run with double precision
    global T = Float64
    global P2 = Point{2,T}
    global P3 = Point{3,T}
    include(testfile)
  end
end
