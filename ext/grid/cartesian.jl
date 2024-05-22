# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function vizgrid2D!(plot::Viz{<:Tuple{CartesianGrid}})
  grid = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colormap = plot[:colormap]
  colorrange = plot[:colorrange]
  showsegments = plot[:showsegments]
  segmentcolor = plot[:segmentcolor]
  segmentsize = plot[:segmentsize]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colormap, $colorrange, $alpha)

  # number of vertices and colors
  nv = Makie.@lift nvertices($grid)
  nc = Makie.@lift $colorant isa AbstractVector ? length($colorant) : 1

  # origin, spacing and size of grid
  or = Makie.@lift ustrip.(to(minimum($grid)))
  sp = Makie.@lift ustrip.(spacing($grid))
  sz = Makie.@lift size($grid)

  if nc[] == 1
    # visualize bounding box with a single
    # color for maximum performance
    bbox = Makie.@lift boundingbox($grid)
    viz!(plot, bbox, color=colorant)

    if showsegments[]
      xyz = Makie.@lift map(x -> ustrip.(x), Meshes.xyz($grid))
      tup = Makie.@lift xysegments($xyz...)
      x, y = Makie.@lift($tup[1]), Makie.@lift($tup[2])
      Makie.lines!(plot, x, y, color=segmentcolor, linewidth=segmentsize)
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

    if showsegments[]
      tup = Makie.@lift xysegments(0:$sz[1], 0:$sz[2])
      x, y = Makie.@lift($tup[1]), Makie.@lift($tup[2])
      Makie.lines!(plot, x, y, color=segmentcolor, linewidth=segmentsize)
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
  colormap = plot[:colormap]
  colorrange = plot[:colorrange]
  showsegments = plot[:showsegments]
  segmentcolor = plot[:segmentcolor]
  segmentsize = plot[:segmentsize]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colormap, $colorrange, $alpha)

  # number of vertices and colors
  nv = Makie.@lift nvertices($grid)
  nc = Makie.@lift $colorant isa AbstractVector ? length($colorant) : 1

  # spacing and coordinates
  sp = Makie.@lift ustrip.(spacing($grid))
  xyz = Makie.@lift map(x -> ustrip.(x), Meshes.xyz($grid))

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

  if showsegments[]
    tup = Makie.@lift xyzsegments($xyz...)
    x, y, z = Makie.@lift($tup[1]), Makie.@lift($tup[2]), Makie.@lift($tup[3])
    Makie.lines!(plot, x, y, z, color=segmentcolor, linewidth=segmentsize)
  end
end
