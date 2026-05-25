# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function vizgrid!(plot::Viz{<:Tuple{CartesianGrid}}, ::Type{<:𝔼}, ::Val{2}, ::Val{2})
  grid = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colormap = plot[:colormap]
  colorrange = plot[:colorrange]
  showsegments = plot[:showsegments]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colormap, $colorrange, $alpha)

  # number of vertices and colors
  nv = Makie.@lift nvertices($grid)
  nc = Makie.@lift $colorant isa AbstractVector ? length($colorant) : 1

  # size and extrema coordinates
  sz = Makie.@lift size($grid)
  xy = Makie.@lift let
    x, y = Meshes.xyz($grid)
    xₛ, xₑ = extrema(ustrip.(x))
    yₛ, yₑ = extrema(ustrip.(y))
    (xₛ, xₑ), (yₛ, yₑ)
  end
  x = Makie.@lift $xy[1]
  y = Makie.@lift $xy[2]

  if nc[] == nv[]
    # visualize as built-in image with interpolation
    C = Makie.@lift reshape($colorant, $sz .+ 1)
    Makie.image!(plot, x, y, C, interpolate=true)
  else
    # visualize as built-in image without interpolation
    C = Makie.@lift $nc == 1 ? fill($colorant, $sz) : reshape($colorant, $sz)
    Makie.image!(plot, x, y, C, interpolate=false)
  end

  if showsegments[]
    vizfacets!(plot)
  end
end

function vizgrid!(plot::Viz{<:Tuple{CartesianGrid}}, ::Type{<:𝔼}, ::Val{3}, ::Val{3})
  # retrieve parameters
  grid = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colormap = plot[:colormap]
  colorrange = plot[:colorrange]
  showsegments = plot[:showsegments]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colormap, $colorrange, $alpha)

  # number of vertices and colors
  nv = Makie.@lift nvertices($grid)
  nc = Makie.@lift $colorant isa AbstractVector ? length($colorant) : 1

  # spacing and coordinates
  sp = Makie.@lift ustrip.(spacing($grid))
  xyz = Makie.@lift map(x -> ustrip.(x), Meshes.xyz($grid))

  if nc[] == nv[]
    # visualize as a quadrangle mesh so that
    # colors can be specified at vertices
    verts = Makie.@lift map(asmakie, eachvertex($grid))
    quads = Makie.@lift reduce(
      vcat,
      map(elements(topology($grid))) do elem
        i1, i2, i3, i4, i5, i6, i7, i8 = indices(elem)
        [
          # follow the same order of vertices specified
          # in the `boundary` method for Hexahedron
          GB.QuadFace(i4, i3, i2, i1),
          GB.QuadFace(i6, i5, i1, i2),
          GB.QuadFace(i3, i7, i6, i2),
          GB.QuadFace(i4, i8, i7, i3),
          GB.QuadFace(i1, i5, i8, i4),
          GB.QuadFace(i6, i7, i8, i5)
        ]
      end
    )
    mesh = Makie.@lift GB.Mesh($verts, $quads)
    Makie.mesh!(plot, mesh, color=colorant, shading=true)
  else
    # visualize as built-in meshscatter
    xs = Makie.@lift $xyz[1][(begin + 1):end]
    ys = Makie.@lift $xyz[2][(begin + 1):end]
    zs = Makie.@lift $xyz[3][(begin + 1):end]
    rect = Makie.@lift Makie.Rect3(-1 .* $sp, $sp)
    coords = Makie.@lift [(x, y, z) for z in $zs for y in $ys for x in $xs]
    Makie.meshscatter!(plot, coords, marker=rect, markersize=1, color=colorant)
  end

  if showsegments[]
    vizfacets!(plot)
  end
end

function vizgridfacets!(plot::Viz{<:Tuple{CartesianGrid}}, ::Type{<:𝔼}, ::Val{2}, ::Val{2})
  grid = plot[:object]
  segmentcolor = plot[:segmentcolor]
  segmentsize = plot[:segmentsize]

  xy = Makie.@lift let
    x, y = Meshes.xyz($grid)
    xysegments(ustrip.(x), ustrip.(y))
  end
  x = Makie.@lift $xy[1]
  y = Makie.@lift $xy[2]

  Makie.lines!(plot, x, y, color=segmentcolor, linewidth=segmentsize)
end

function vizgridfacets!(plot::Viz{<:Tuple{CartesianGrid}}, ::Type{<:𝔼}, ::Val{3}, ::Val{3})
  grid = plot[:object]
  segmentcolor = plot[:segmentcolor]
  segmentsize = plot[:segmentsize]

  xyz = Makie.@lift let
    x, y, z = Meshes.xyz($grid)
    xyzsegments(ustrip.(x), ustrip.(y), ustrip.(z))
  end
  x = Makie.@lift $xyz[1]
  y = Makie.@lift $xyz[2]
  z = Makie.@lift $xyz[3]

  Makie.lines!(plot, x, y, z, color=segmentcolor, linewidth=segmentsize)
end
