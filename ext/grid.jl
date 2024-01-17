# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{Grid}})
  grid = plot[:object][]
  T = typeof(grid)
  Dim = embeddim(grid)
  if Dim == 1
    vizgrid1D!(T, plot)
  elseif Dim == 2
    vizgrid2D!(T, plot)
  elseif Dim == 3
    vizgrid3D!(T, plot)
  end
end

vizgrid1D!(::Type{<:Grid}, plot) = vizmesh1D!(plot)

vizgrid2D!(::Type{<:Grid}, plot) = vizmesh2D!(plot)

function vizgrid2D!(::Type{<:Union{CartesianGrid,RectilinearGrid}}, plot)
  grid = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colorscheme = plot[:colorscheme]
  segmentsize = plot[:segmentsize]
  showfacets = plot[:showfacets]
  facetcolor = plot[:facetcolor]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colorscheme, $alpha)

  # number of vertices and colors
  nv = Makie.@lift nvertices($grid)
  nc = Makie.@lift $colorant isa AbstractVector ? length($colorant) : 1

  # grid coordinates
  xyz = Makie.@lift Meshes.xyz($grid)
  xs = Makie.@lift $xyz[1]
  ys = Makie.@lift $xyz[2]

  if nc[] == 1
    # visualize bounding box with a single
    # color for maximum performance
    bbox = Makie.@lift boundingbox($grid)
    viz!(plot, bbox, color=colorant)
  else
    if nc[] == nv[]
      # visualize as a simple mesh so that
      # colors can be specified at vertices
      vizmesh2D!(plot)
    else
      # visualize as built-in heatmap
      sz = Makie.@lift size($grid)
      C = Makie.@lift reshape($colorant, $sz)
      Makie.heatmap!(plot, xs, ys, C)
    end
  end

  if showfacets[]
    tup = Makie.@lift xysegments($xs, $ys)
    x, y = Makie.@lift($tup[1]), Makie.@lift($tup[2])
    Makie.lines!(plot, x, y, color=facetcolor, linewidth=segmentsize)
  end
end

function vizgrid3D!(::Type{<:Grid}, plot)
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

function vizgrid3D!(::Type{<:CartesianGrid}, plot)
  # retrieve parameters
  grid = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colorscheme = plot[:colorscheme]
  segmentsize = plot[:segmentsize]
  showfacets = plot[:showfacets]
  facetcolor = plot[:facetcolor]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colorscheme, $alpha)

  # number of vertices and colors
  nv = Makie.@lift nvertices($grid)
  nc = Makie.@lift $colorant isa AbstractVector ? length($colorant) : 1

  # spacing and coordinates of grid
  sp = Makie.@lift spacing($grid)
  xyz = Makie.@lift Meshes.xyz($grid)

  if nc[] == 1
    # visualize bounding box with a single
    # color for maximum performance
    bbox = Makie.@lift boundingbox($grid)
    viz!(plot, bbox, color=colorant)
  else
    if nc[] == nv[]
      error("not implemented")
    else
      # visualize as built-in meshscatter
      xs = Makie.@lift $xyz[1][(begin + 1):end]
      ys = Makie.@lift $xyz[2][(begin + 1):end]
      zs = Makie.@lift $xyz[3][(begin + 1):end]
      rect = Makie.@lift Makie.Rect3(-1 .* $sp, $sp)
      coords = Makie.@lift [(x, y, z) for z in $zs for y in $ys for x in $xs]
      Makie.meshscatter!(plot, coords, marker=rect, markersize=1, color=colorant)
    end
  end

  if showfacets[]
    tup = Makie.@lift xyzsegments($xyz...)
    x, y, z = Makie.@lift($tup[1]), Makie.@lift($tup[2]), Makie.@lift($tup[3])
    Makie.lines!(plot, x, y, z, color=facetcolor, linewidth=segmentsize)
  end
end

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
