using TestItems
using TestItemRunner

@run_package_tests

@testsnippet Setup begin
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
  using StableRNGs
  using ReferenceTests, ImageIO

  using TransformsBase: Identity, →

  import TransformsBase as TB
  import CairoMakie as Mke

  # environment settings
  isCI = "CI" ∈ keys(ENV)
  islinux = Sys.islinux()
  visualtests = !isCI || (isCI && islinux)
  datadir = joinpath(@__DIR__, "data")

  # float settings
  T = Float64

  include("testutils.jl")
end
