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
  using TransformsBase
  using DelimitedFiles
  using Unitful
  using Rotations
  using StableRNGs

  import TransformsBase as TB
  import Mooncake

  # environment settings
  isCI = "CI" ∈ keys(ENV)
  datadir = joinpath(@__DIR__, "data")

  # float settings
  T = if isCI
    if ENV["FLOAT_TYPE"] == "Float32"
      Float32
    elseif ENV["FLOAT_TYPE"] == "Float64"
      Float64
    end
  else
    Float64
  end

  include("testutils.jl")
end
