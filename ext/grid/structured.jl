# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function vizgrid!(plot::Viz{<:Tuple{StructuredGrid}}, M::Type{<:𝔼}, pdim::Val{2}, edim::Val{2})
  showsegments = plot.showsegments

  if crs(plot.object[]) <: Cartesian
    # process color spec into colorant
    Makie.map!(process, plot, [:color, :colormap, :colorrange, :alpha], :colorant)

    # number of vertices and colors
    Makie.map!(nvertices, plot, [:object], :nv)
    Makie.map!(plot, [:colorant], :nc) do colorant
      colorant isa AbstractVector ? length(colorant) : 1
    end

    if plot.nc[] == plot.nv[]
      # size and coordinates
      Makie.map!(plot, [:object], :sz) do grid
        size(grid) .+ 1
      end
      Makie.map!(plot, [:object], [:X, :Y]) do grid
        map(X -> ustrip.(X), Meshes.XYZ(grid))
      end

      # visualize as built-in surface
      Makie.map!(plot, [:colorant, :sz], :C) do colorant, sz
        reshape(colorant, sz)
      end
      Makie.surface!(plot, plot.X, plot.Y, color=plot.C)

      if showsegments[]
        vizfacets!(plot)
      end
    else
      vizmesh!(plot)
    end
  else
    vizgridfallback!(plot, M, pdim, edim)
  end
end

function vizgridfacets!(plot::Viz{<:Tuple{StructuredGrid}}, ::Type{<:𝔼}, ::Val{2}, ::Val{2})
  segmentcolor = plot.segmentcolor
  segmentsize = plot.segmentsize

  Makie.map!(structuredsegments, plot, [:object], :facets_tup)
  Makie.map!(plot, [:facets_tup], [:facets_x, :facets_y]) do tup
    (tup[1], tup[2])
  end
  Makie.lines!(plot, plot.facets_x, plot.facets_y, color=segmentcolor, linewidth=segmentsize)
end

function structuredsegments(grid)
  cinds = CartesianIndices(size(grid) .+ 1)
  coords = []
  # vertical segments
  for j in axes(cinds, 2)
    for i in axes(cinds, 1)
      p = vertex(grid, cinds[i, j])
      push!(coords, Tuple(ustrip.(to(p))))
    end
    push!(coords, (NaN, NaN))
  end
  # horizontal segments
  for i in axes(cinds, 1)
    for j in axes(cinds, 2)
      p = vertex(grid, cinds[i, j])
      push!(coords, Tuple(ustrip.(to(p))))
    end
    push!(coords, (NaN, NaN))
  end
  x = getindex.(coords, 1)
  y = getindex.(coords, 2)
  x, y
end
