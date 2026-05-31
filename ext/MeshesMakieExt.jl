# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

module MeshesMakieExt

using Meshes
using Unitful
using Rotations
using StaticArrays
using LinearAlgebra
using CoordRefSystems
using Colorfy

using Unitful: numtype
using Meshes: lentype

import TransformsBase as TB
import Makie.GeometryBasics as GB

import Meshes: viz, viz!
import Makie

Makie.@recipe(Viz, object) do scene
  Makie.Attributes(
    color=:slategray3,
    alpha=nothing,
    colormap=nothing,
    colorrange=nothing,
    showsegments=false,
    segmentcolor=:gray30,
    segmentsize=1.5,
    showpoints=false,
    pointmarker=:circle,
    pointcolor=:gray30,
    pointsize=4
  )
end

# choose between 2D and 3D axis
Makie.args_preferred_axis(g::Geometry) = axis(embeddim(g))
Makie.args_preferred_axis(d::Domain) = axis(embeddim(d))
Makie.args_preferred_axis(::Vec{Dim}) where {Dim} = axis(Dim)
Makie.args_preferred_axis(g::AbstractVector{<:Geometry}) = Makie.args_preferred_axis(first(g))
Makie.args_preferred_axis(v::AbstractVector{<:Vec}) = Makie.args_preferred_axis(first(v))

axis(dim) = dim === 3 ? Makie.Axis3 : Makie.Axis

# choose default axis attributes
Makie.preferred_axis_attributes(_, g::Geometry) = axisattributes(g)
Makie.preferred_axis_attributes(_, d::Domain) = axisattributes(d)
Makie.preferred_axis_attributes(_, ::Vec{Dim,ℒ}) where {Dim,ℒ} = axisattributes(Dim, ℒ)
Makie.preferred_axis_attributes(A, g::AbstractVector{<:Geometry}) = Makie.preferred_axis_attributes(A, first(g))
Makie.preferred_axis_attributes(A, v::AbstractVector{<:Vec}) = Makie.preferred_axis_attributes(A, first(v))

function axisattributes(g)
  aspect = fixaspect(embeddim(g))
  labels = xyzlabels(crs(g))
  merge(aspect, labels)
end

function axisattributes(Dim, ℒ)
  aspect = fixaspect(Dim)
  labels = xyzlabels(Cartesian{NoDatum,Dim,ℒ})
  merge(aspect, labels)
end

fixaspect(dim) = dim === 3 ? (aspect=:data, viewmode=:free, perspectiveness=0.5) : (aspect=Makie.DataAspect(),)

function xyzlabels(CRS)
  u = CoordRefSystems.lentype(CRS) |> unit
  if CoordRefSystems.ndims(CRS) === 3
    (xlabel="X [$u]", ylabel="Y [$u]", zlabel="Z [$u]")
  elseif CRS <: CoordRefSystems.Projected
    (xlabel="Easting [$u]", ylabel="Northing [$u]")
  else
    (xlabel="X [$u]", ylabel="Y [$u]")
  end
end

# color handling
include("colors.jl")

# utilities
include("utils.jl")

# viz recipes
include("mesh.jl")
include("grid.jl")
include("geomset.jl")
include("subdomain.jl")
include("fallbacks.jl")

end
