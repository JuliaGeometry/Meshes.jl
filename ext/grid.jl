# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{Grid}})
  grid = plot[:object][]
  Dim = embeddim(grid)
  if Dim == 2
    vizgrid2D!(plot)
  elseif Dim == 3
    vizgrid3D!(plot)
  end
end

vizgrid2D!(plot) = vizmesh2D!(plot)

function vizgrid3D!(plot)
  grid = plot[:object]
  color = plot[:color]

  # number of vertices and colors
  nv = Makie.@lift nvertices($grid)
  nc = Makie.@lift $color isa AbstractVector ? length($color) : 1

  if nv[] == nc[]
    error("not implemented")
  else
    vizmesh3D!(plot)
  end
end

# defining a Makie.data_limits method is necessary because
# Makie.scale!, Makie.translate! and Makie.rotate!
# don't adjust axis limits automatically
function Makie.data_limits(plot::Viz{<:Tuple{Grid}})
  grid = plot[:object][]
  bbox = boundingbox(grid)
  pmin = aspoint3f(minimum(bbox))
  pmax = aspoint3f(maximum(bbox))
  Makie.Rect3f([pmin, pmax])
end

aspoint3f(p::Point{2}) = Makie.Point3f(ustrip.(coordinates(p))..., 0)
aspoint3f(p::Point{3}) = Makie.Point3f(ustrip.(coordinates(p))...)

# ----------------
# SPECIALIZATIONS
# ----------------

include("grid/cartesian.jl")
include("grid/rectilinear.jl")
include("grid/transformed.jl")
include("grid/structured.jl")

# helper functions to create a minimum number
# of line segments within Cartesian/Rectilinear grid
function xysegments(xs, ys)
  coords = []
  for x in xs
    push!(coords, (x, first(ys)))
    push!(coords, (x, last(ys)))
    push!(coords, (NaN, NaN))
  end
  for y in ys
    push!(coords, (first(xs), y))
    push!(coords, (last(xs), y))
    push!(coords, (NaN, NaN))
  end
  x = getindex.(coords, 1)
  y = getindex.(coords, 2)
  x, y
end

function xyzsegments(xs, ys, zs)
  coords = []
  for y in ys, z in zs
    push!(coords, (first(xs), y, z))
    push!(coords, (last(xs), y, z))
    push!(coords, (NaN, NaN, NaN))
  end
  for x in xs, z in zs
    push!(coords, (x, first(ys), z))
    push!(coords, (x, last(ys), z))
    push!(coords, (NaN, NaN, NaN))
  end
  for x in xs, y in ys
    push!(coords, (x, y, first(zs)))
    push!(coords, (x, y, last(zs)))
    push!(coords, (NaN, NaN, NaN))
  end
  x = getindex.(coords, 1)
  y = getindex.(coords, 2)
  z = getindex.(coords, 3)
  x, y, z
end
