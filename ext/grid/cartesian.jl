# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function vizgrid2D!(plot::Viz{<:Tuple{CartesianGrid}})
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

  # origin, spacing and size of grid
  or = Makie.@lift coordinates(minimum($grid))
  sp = Makie.@lift spacing($grid)
  sz = Makie.@lift size($grid)

  if nc[] == 1
    # visualize bounding box with a single
    # color for maximum performance
    bbox = Makie.@lift boundingbox($grid)
    viz!(plot, bbox, color=colorant)

    if showfacets[]
      tup = Makie.@lift xysegments(Meshes.xyz($grid)...)
      x, y = Makie.@lift($tup[1]), Makie.@lift($tup[2])
      Makie.lines!(plot, x, y, color=facetcolor, linewidth=segmentsize)
    end
  else
    if nc[] == nv[]
      # visualize as built-in image with interpolation
      C = Makie.@lift reshape($colorant, $sz .+ 1)
      Makie.image!(plot, C, interpolate=true)
    else
      # visualize as built-in image without interpolation
      C = Makie.@lift reshape($colorant, $sz)
      Makie.image!(plot, C, interpolate=false)
    end

    if showfacets[]
      tup = Makie.@lift xysegments(0:$sz[1], 0:$sz[2])
      x, y = Makie.@lift($tup[1]), Makie.@lift($tup[2])
      Makie.lines!(plot, x, y, color=facetcolor, linewidth=segmentsize)
    end

    # adjust spacing and origin
    spx, spy = sp[]
    orx, ory = or[]
    Makie.scale!(plot, spx, spy)
    Makie.translate!(plot, orx, ory)
  end
end

function vizgrid3D!(plot::Viz{<:Tuple{CartesianGrid}})
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

  # spacing and coordinates
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
