# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function vizgrid!(plot::Viz{<:Tuple{CartesianGrid}}, ::Type{<:𝔼}, ::Val{2}, ::Val{2})
  showsegments = plot[:showsegments]

  # process color spec into colorant
  Makie.map!(process, plot, [:color, :colormap, :colorrange, :alpha], :colorant)

  # number of vertices and colors
  Makie.map!(nvertices, plot, [:object], :nv)
  Makie.map!(plot, [:colorant], :nc) do colorant
    colorant isa AbstractVector ? length(colorant) : 1
  end

  # size and extrema coordinates
  Makie.map!(size, plot, [:object], :sz)
  Makie.map!(plot, [:object], :xy) do grid
    x, y = Meshes.xyz(grid)
    xₛ, xₑ = extrema(ustrip.(x))
    yₛ, yₑ = extrema(ustrip.(y))
    (xₛ, xₑ), (yₛ, yₑ)
  end
  Makie.map!(plot, [:xy], :x) do xy
    xy[1]
  end
  Makie.map!(plot, [:xy], :y) do xy
    xy[2]
  end

  if plot[:nc][] == plot[:nv][]
    # visualize as built-in image with interpolation
    Makie.map!(plot, [:colorant, :sz], :C) do colorant, sz
      reshape(colorant, sz .+ 1)
    end
    Makie.image!(plot, plot[:x], plot[:y], plot[:C], interpolate=true)
  else
    # visualize as built-in image without interpolation
    Makie.map!(plot, [:nc, :colorant, :sz], :C) do nc, colorant, sz
      nc == 1 ? fill(colorant, sz) : reshape(colorant, sz)
    end
    Makie.image!(plot, plot[:x], plot[:y], plot[:C], interpolate=false)
  end

  if showsegments[]
    vizfacets!(plot)
  end
end

function vizgrid!(plot::Viz{<:Tuple{CartesianGrid}}, ::Type{<:𝔼}, ::Val{3}, ::Val{3})
  # retrieve parameters
  showsegments = plot[:showsegments]

  # process color spec into colorant
  Makie.map!(process, plot, [:color, :colormap, :colorrange, :alpha], :colorant)

  # number of vertices and colors
  Makie.map!(nvertices, plot, [:object], :nv)
  Makie.map!(plot, [:colorant], :nc) do colorant
    colorant isa AbstractVector ? length(colorant) : 1
  end

  # spacing and coordinates
  Makie.map!(plot, [:object], :sp) do grid
    ustrip.(spacing(grid))
  end
  Makie.map!(plot, [:object], :xyz) do grid
    map(x -> ustrip.(x), Meshes.xyz(grid))
  end

  if plot[:nc][] == plot[:nv][]
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
    Makie.mesh!(plot, plot[:mesh], color=plot[:colorant], shading=true)
  else
    # visualize as built-in meshscatter
    Makie.map!(plot, [:xyz], :xs) do xyz
      xyz[1][(begin + 1):end]
    end
    Makie.map!(plot, [:xyz], :ys) do xyz
      xyz[2][(begin + 1):end]
    end
    Makie.map!(plot, [:xyz], :zs) do xyz
      xyz[3][(begin + 1):end]
    end
    Makie.map!(plot, [:sp], :rect) do sp
      Makie.Rect3(-1 .* sp, sp)
    end
    Makie.map!(plot, [:xs, :ys, :zs], :coords) do xs, ys, zs
      [(x, y, z) for z in zs for y in ys for x in xs]
    end
    Makie.meshscatter!(plot, plot[:coords], marker=plot[:rect], markersize=1, color=plot[:colorant])
  end

  if showsegments[]
    vizfacets!(plot)
  end
end

function vizgridfacets!(plot::Viz{<:Tuple{CartesianGrid}}, ::Type{<:𝔼}, ::Val{2}, ::Val{2})
  segmentcolor = plot[:segmentcolor]
  segmentsize = plot[:segmentsize]

  Makie.map!(plot, [:object], :xy) do grid
    x, y = Meshes.xyz(grid)
    xysegments(ustrip.(x), ustrip.(y))
  end
  Makie.map!(plot, [:xy], :x) do xy
    xy[1]
  end
  Makie.map!(plot, [:xy], :y) do xy
    xy[2]
  end

  Makie.lines!(plot, plot[:x], plot[:y], color=segmentcolor, linewidth=segmentsize)
end

function vizgridfacets!(plot::Viz{<:Tuple{CartesianGrid}}, ::Type{<:𝔼}, ::Val{3}, ::Val{3})
  segmentcolor = plot[:segmentcolor]
  segmentsize = plot[:segmentsize]

  Makie.map!(plot, [:object], :xyz) do grid
    x, y, z = Meshes.xyz(grid)
    xyzsegments(ustrip.(x), ustrip.(y), ustrip.(z))
  end
  Makie.map!(plot, [:xyz], :x) do xyz
    xyz[1]
  end
  Makie.map!(plot, [:xyz], :y) do xyz
    xyz[2]
  end
  Makie.map!(plot, [:xyz], :z) do xyz
    xyz[3]
  end

  Makie.lines!(plot, plot[:x], plot[:y], plot[:z], color=segmentcolor, linewidth=segmentsize)
end
