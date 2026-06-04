# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function vizgrid!(plot::Viz{<:Tuple{CartesianGrid}}, ::Type{<:𝔼}, ::Val{2}, ::Val{2})
  showsegments = plot.showsegments

  # process color spec into colorant
  Makie.map!(process, plot, [:color, :colormap, :colorrange, :alpha], :colorant)

  # number of vertices and colors
  Makie.map!(nvertices, plot, [:object], :nv)
  Makie.map!(plot, [:colorant], :nc) do colorant
    colorant isa AbstractVector ? length(colorant) : 1
  end

  # size and extrema coordinates
  Makie.map!(size, plot, [:object], :sz)
  Makie.map!(plot, [:object], [:x, :y]) do grid
    x, y = Meshes.xyz(grid)
    xₛ, xₑ = extrema(ustrip.(x))
    yₛ, yₑ = extrema(ustrip.(y))
    ((xₛ, xₑ), (yₛ, yₑ))
  end

  if plot.nc[] == plot.nv[]
    # visualize as built-in image with interpolation
    Makie.map!(plot, [:colorant, :sz], :C) do colorant, sz
      reshape(colorant, sz .+ 1)
    end
    Makie.image!(plot, plot.x, plot.y, plot.C, interpolate=true)
  else
    # visualize as built-in image without interpolation
    Makie.map!(plot, [:colorant, :sz, :nc], :C) do colorant, sz, nc
      nc == 1 ? fill(colorant, sz) : reshape(colorant, sz)
    end
    Makie.image!(plot, plot.x, plot.y, plot.C, interpolate=false)
  end

  if showsegments[]
    vizfacets!(plot)
  end
end

function vizgrid!(plot::Viz{<:Tuple{CartesianGrid}}, ::Type{<:𝔼}, ::Val{3}, ::Val{3})
  # retrieve parameters
  showsegments = plot.showsegments

  # process color spec into colorant
  Makie.map!(process, plot, [:color, :colormap, :colorrange, :alpha], :colorant)

  # number of vertices and colors
  Makie.map!(nvertices, plot, [:object], :nv)
  Makie.map!(plot, [:colorant], :nc) do colorant
    colorant isa AbstractVector ? length(colorant) : 1
  end

  if plot.nc[] == plot.nv[]
    # visualize as a quadrangle mesh so that
    # colors can be specified at vertices
    Makie.map!(plot, [:object], :verts) do grid
      map(asmakie, eachvertex(grid))
    end
    Makie.map!(plot, [:object], :quads) do grid
      reduce(
        vcat,
        map(elements(topology(grid))) do elem
          i1, i2, i3, i4, i5, i6, i7, i8 = indices(elem)
          [
            # follow the same order of vertices specified
            # in the `boundary` method for `Hexahedron`
            GB.QuadFace(i4, i3, i2, i1),
            GB.QuadFace(i6, i5, i1, i2),
            GB.QuadFace(i3, i7, i6, i2),
            GB.QuadFace(i4, i8, i7, i3),
            GB.QuadFace(i1, i5, i8, i4),
            GB.QuadFace(i6, i7, i8, i5)
          ]
        end
      )
    end
    Makie.map!(GB.Mesh, plot, [:verts, :quads], :mesh)
    Makie.mesh!(plot, plot.mesh, color=plot.colorant, shading=true)
  else
    # visualize as built-in meshscatter
    Makie.map!(plot, [:object], [:coords, :rect]) do grid
      xyz = map(x -> ustrip.(x), Meshes.xyz(grid))
      xs = xyz[1][(begin + 1):end]
      ys = xyz[2][(begin + 1):end]
      zs = xyz[3][(begin + 1):end]
      coords = [(x, y, z) for z in zs for y in ys for x in xs]
      sp = ustrip.(spacing(grid))
      rect = Makie.Rect3(-1 .* sp, sp)
      (coords, rect)
    end
    Makie.meshscatter!(plot, plot.coords, marker=plot.rect, markersize=1, color=plot.colorant)
  end

  if showsegments[]
    vizfacets!(plot)
  end
end

function vizgridfacets!(plot::Viz{<:Tuple{CartesianGrid}}, ::Type{<:𝔼}, ::Val{2}, ::Val{2})
  segmentcolor = plot.segmentcolor
  segmentsize = plot.segmentsize

  Makie.map!(plot, [:object], [:facets_x, :facets_y]) do grid
    x, y = Meshes.xyz(grid)
    xysegments(ustrip.(x), ustrip.(y))
  end

  Makie.lines!(plot, plot.facets_x, plot.facets_y, color=segmentcolor, linewidth=segmentsize)
end

function vizgridfacets!(plot::Viz{<:Tuple{CartesianGrid}}, ::Type{<:𝔼}, ::Val{3}, ::Val{3})
  segmentcolor = plot.segmentcolor
  segmentsize = plot.segmentsize

  Makie.map!(plot, [:object], [:facets_x, :facets_y, :facets_z]) do grid
    x, y, z = Meshes.xyz(grid)
    xyzsegments(ustrip.(x), ustrip.(y), ustrip.(z))
  end

  Makie.lines!(plot, plot.facets_x, plot.facets_y, plot.facets_z, color=segmentcolor, linewidth=segmentsize)
end
