# Setup code has to be in this quote
# so it can be executed twice (for f32 and f64)

run_setup = quote

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
  using Rotations
  using Test, Random
  using ReferenceTests, ImageIO

  import TransformsBase as TB
  import CairoMakie as Mke

  # environment settings
  isCI = "CI" ∈ keys(ENV)
  islinux = Sys.islinux()
  visualtests = !isCI || (isCI && islinux)
  datadir = joinpath(@__DIR__, "data")

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
      @assert first(outer) == last(outer)
      @assert all(first(i) == last(i) for i in inners)
      PolyArea(outer[begin:(end - 1)], [i[begin:(end - 1)] for i in inners])
    end
  end

  # helper function to read *.ply files containing meshes
  function readply(T, fname)
    ply = load_ply(fname)
    x = ply["vertex"]["x"]
    y = ply["vertex"]["y"]
    z = ply["vertex"]["z"]
    points = Point{3,T}.(x, y, z)
    connec = [connect(Tuple(c .+ 1)) for c in ply["face"]["vertex_indices"]]
    SimpleMesh(points, connec)
  end

  include("dummy.jl")
  include("allocation_macros.jl")

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
    "sets.jl",
    "mesh.jl",
    "data.jl",
    "meshdata.jl",
    "utils.jl",
    "views.jl",
    "viewing.jl",
    "partitioning.jl",
    "traversing.jl",
    "neighborhoods.jl",
    "neighborsearch.jl",
    "distances.jl",
    "supportfun.jl",
    "matrices.jl",
    "merging.jl",
    "sampling.jl",
    "intersections.jl",
    "pointification.jl",
    "discretization.jl",
    "simplification.jl",
    "refinement.jl",
    "boundingboxes.jl",
    "hulls.jl",
    "transforms.jl"
  ]

end # end quote


# --------------------------------
# RUN TESTS WITH SINGLE PRECISION
# --------------------------------
module TestsSinglePrecision
  using ..Main: run_setup

  eval(run_setup)

  const T::Type = Float32
  const P1::Type, P2::Type, P3::Type = Point{1,T}, Point{2,T}, Point{3,T}
  const V1::Type, V2::Type, V3::Type = Vec{1,T}, Vec{2,T}, Vec{3,T}
  @testset "Meshes.jl ($T)" begin
    for testfile in testfiles
      println("Testing $testfile...")
      include(testfile)
    end
  end

end


# --------------------------------
# RUN TESTS WITH DOUBLE PRECISION
# --------------------------------
module TestsDoublePrecision
  using ..Main: run_setup

  eval(run_setup)

  const T::Type = Float64
  const P1::Type, P2::Type, P3::Type = Point{1,T}, Point{2,T}, Point{3,T}
  const V1::Type, V2::Type, V3::Type = Vec{1,T}, Vec{2,T}, Vec{3,T}
  @testset "Meshes.jl ($T)" begin
    for testfile in testfiles
      println("Testing $testfile...")
      include(testfile)
    end
  end
end
