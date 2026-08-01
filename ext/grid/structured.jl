# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function vizgrid!(plot::Viz{<:Tuple{StructuredGrid}}, M::Type{<:𝔼}, pdim::Val{2}, edim::Val{2})
  if crs(plot.object[]) <: Cartesian
    # number of vertices and colors
    Makie.map!(plot, [:object, :colorant], [:nv, :nc]) do grid, colorant
      nv = nvertices(grid)
      nc = colorant isa AbstractVector ? length(colorant) : 1
      nv, nc
    end

    if plot.nc[] == plot.nv[]
      # visualize as built-in surface
      Makie.map!(plot, [:object, :colorant], [:X, :Y, :C, :usetiles]) do grid, colorant
        ne = nelements(grid)
        X, Y = map(c -> ustrip.(c), Meshes.XYZ(grid))
        C = reshape(colorant, Meshes.vsize(grid))
        X, Y, C, (ne ≥ 100000)
      end

      if plot.usetiles[]
        for (i, tile) in enumerate(TileIterator(axes(plot.C[]), (100, 100)))
          Xi = Symbol(:X, i)
          Yi = Symbol(:Y, i)
          Ci = Symbol(:C, i)
          Makie.map!(plot, [:X, :Y, :C], [Xi, Yi, Ci]) do X, Y, C
            view(X, tile...), view(Y, tile...), view(C, tile...)
          end
          Makie.surface!(plot, plot[Xi], plot[Yi], color=plot[Ci])
        end
      else
        Makie.surface!(plot, plot.X, plot.Y, color=plot.C)
      end

      if plot.showsegments[]
        vizfacets!(plot)
      end
    else
      # visualize as simple mesh
      vizmesh!(plot)
    end
  else
    vizgridfallback!(plot, M, pdim, edim)
  end
end

function vizgridfacets!(plot::Viz{<:Tuple{StructuredGrid}}, ::Type{<:𝔼}, ::Val{2}, ::Val{2})
  Makie.map!(structuredsegments, plot, :object, [:xfacets, :yfacets])
  Makie.lines!(plot, plot.xfacets, plot.yfacets, color=plot.segmentcolor, linewidth=plot.segmentsize)
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
