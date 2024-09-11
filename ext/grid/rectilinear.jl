# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function vizgrid!(plot::Viz{<:Tuple{RectilinearGrid}}, M::Type{<:𝔼}, pdim::Val{2}, edim::Val{2})
  grid = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colormap = plot[:colormap]
  colorrange = plot[:colorrange]
  showsegments = plot[:showsegments]
  segmentcolor = plot[:segmentcolor]
  segmentsize = plot[:segmentsize]

  if crs(grid[]) <: Cartesian
    # process color spec into colorant
    colorant = Makie.@lift process($color, $colormap, $colorrange, $alpha)

    # number of vertices and colors
    nv = Makie.@lift nvertices($grid)
    nc = Makie.@lift $colorant isa AbstractVector ? length($colorant) : 1

    # grid coordinates
    xyz = Makie.@lift map(x -> ustrip.(x), Meshes.xyz($grid))
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
        vizmesh!(plot, M, pdim, edim)
      else
        # visualize as built-in heatmap
        sz = Makie.@lift size($grid)
        C = Makie.@lift reshape($colorant, $sz)
        Makie.heatmap!(plot, xs, ys, C)
      end
    end

    if showsegments[]
      tup = Makie.@lift xysegments($xs, $ys)
      x, y = Makie.@lift($tup[1]), Makie.@lift($tup[2])
      Makie.lines!(plot, x, y, color=segmentcolor, linewidth=segmentsize)
    end
  else
    vizgridfallback!(plot, M, pdim, edim)
  end
end
