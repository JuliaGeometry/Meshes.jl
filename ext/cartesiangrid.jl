# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

Makie.plottype(::CartesianGrid) = Viz{<:Tuple{CartesianGrid}}

function Makie.plot!(plot::Viz{<:Tuple{CartesianGrid}})
  # retrieve parameters
  grid = plot[:object][]
  color = plot[:color][]
  showfacets = plot[:showfacets][]
  ndim = embeddim(grid)

  # different recipes for Cartesian grids
  # with 1D, 2D, 3D elements
  if color isa AbstractVector
    # visualize grid as heatmap or volume
    if ndim == 1
      vizgrid1D!(plot)
    elseif ndim == 2
      if length(color) == nvertices(grid)
        vizmesh2D!(plot)
      else
        vizgrid2D!(plot)
      end
    elseif ndim == 3
      vizgrid3D!(plot)
    end
  else
    # create the smallest mesh of simplices
    vizgrid!(plot)
  end

  if showfacets
    # create minimum number of segments
    vizsegs!(plot)
  end
end

function vizgrid1D!(plot)
  grid = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colorscheme = plot[:colorscheme]
  showfacets = plot[:showfacets]
  facetcolor = plot[:facetcolor]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colorscheme, $alpha)

  cparams = Makie.@lift let
    nd = embeddim($grid)
    or = coordinates(minimum($grid))
    sp = spacing($grid)
    sz = size($grid)

    xs, ys = cartesiancenters(or, sp, sz, nd)

    xs⁻ = [(xs .- sp[1] / 2); (last(xs) + sp[1] / 2)]
    ys⁻ = [ys; last(ys)]

    points = [Point(x, y) for (x, y) in zip(xs⁻, ys⁻)]
    mesh = SimpleMesh(points, topology($grid))

    colors = [$colorant; last($colorant)]

    mesh, colors
  end

  # unpack observable of parameters
  mesh = Makie.@lift $cparams[1]
  colors = Makie.@lift $cparams[2]

  # rely on recipe for simplices
  viz!(plot, mesh, color=colors, showfacets=showfacets, facetcolor=facetcolor)
end

function vizgrid2D!(plot)
  grid = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colorscheme = plot[:colorscheme]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colorscheme, $alpha)

  cparams = Makie.@lift let
    nd = embeddim($grid)
    or = coordinates(minimum($grid))
    sp = spacing($grid)
    sz = size($grid)

    xs, ys = cartesiancenters(or, sp, sz, nd)

    C = reshape($colorant, sz)

    xs, ys, C
  end

  # unpack observable of parameters
  xs = Makie.@lift $cparams[1]
  ys = Makie.@lift $cparams[2]
  C = Makie.@lift $cparams[3]

  Makie.heatmap!(plot, xs, ys, C)
end

function vizgrid3D!(plot)
  grid = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colorscheme = plot[:colorscheme]

  # process color spec into colorant
  colorant = Makie.@lift process($color, $colorscheme, $alpha)

  cparams = Makie.@lift let
    nd = embeddim($grid)
    or = coordinates(minimum($grid))
    sp = spacing($grid)
    sz = size($grid)

    xs, ys, zs = cartesiancenters(or, sp, sz, nd)

    xs⁺ = xs .+ sp[1] / 2
    ys⁺ = ys .+ sp[2] / 2
    zs⁺ = zs .+ sp[3] / 2

    coords = [(x, y, z) for z in zs⁺ for y in ys⁺ for x in xs⁺]

    marker = Makie.Rect3(-1 .* sp, sp)

    coords, marker
  end

  # unpack observable of parameters
  coords = Makie.@lift $cparams[1]
  marker = Makie.@lift $cparams[2]

  Makie.meshscatter!(plot, coords, marker=marker, markersize=1, color=colorant)
end

function vizgrid!(plot)
  grid = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]

  mesh = Makie.@lift let
    nd = embeddim($grid)
    or = coordinates(minimum($grid))
    sp = spacing($grid)
    sz = size($grid)

    cartesianmesh(or, sp, sz, nd)
  end

  viz!(plot, mesh, color=color, alpha=alpha, showfacets=false)
end

function vizsegs!(plot)
  grid = plot[:object]
  facetcolor = plot[:facetcolor]
  segmentsize = plot[:segmentsize]

  cparams = Makie.@lift let
    nd = embeddim($grid)
    or = coordinates(minimum($grid))
    sp = spacing($grid)
    sz = size($grid)

    cartesiansegments(or, sp, sz, nd)
  end

  # unpack observable of parameters
  xyz = [(Makie.@lift $cparams[i]) for i in 1:embeddim(grid[])]

  Makie.lines!(plot, xyz..., color=facetcolor, linewidth=segmentsize)
end

# helper function to create the smallest mesh
# of simplices covering the Cartesian grid
function cartesianmesh(or, sp, sz, nd)
  if nd == 1
    A = Point2(or[1], 0) + Vec2(0, 0)
    B = Point2(or[1], 0) + Vec2(sp[1] * sz[1], 0)
    points = [A, B]
    connec = connect.([(1, 2)])
    SimpleMesh(points, connec)
  elseif nd == 2
    A = Point2(or) + Vec2(0, 0)
    B = Point2(or) + Vec2(sp[1] * sz[1], sp[2] * sz[2])
    discretize(Box(A, B))
  elseif nd == 3
    A = Point3(or) + Vec3(0, 0, 0)
    B = Point3(or) + Vec3(sp[1] * sz[1], sp[2] * sz[2], sp[3] * sz[3])
    boundary(Box(A, B))
  else
    throw(ErrorException("not implemented"))
  end
end

# helper function to create a minimum number
# of line segments within Cartesian grid
function cartesiansegments(or, sp, sz, nd)
  if nd == 2
    xs = range(or[1], step=sp[1], length=sz[1] + 1)
    ys = range(or[2], step=sp[2], length=sz[2] + 1)
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
  elseif nd == 3
    xs = range(or[1], step=sp[1], length=sz[1] + 1)
    ys = range(or[2], step=sp[2], length=sz[2] + 1)
    zs = range(or[3], step=sp[3], length=sz[3] + 1)
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
  else
    throw(ErrorException("not implemented"))
  end
end

# helper function to create the center of
# the elements of the Cartesian grid
function cartesiancenters(or, sp, sz, nd)
  if nd == 1
    xs = range(or[1] + sp[1] / 2, step=sp[1], length=sz[1])
    ys = fill(0.0, sz[1])
    xs, ys
  elseif nd == 2
    xs = range(or[1] + sp[1] / 2, step=sp[1], length=sz[1])
    ys = range(or[2] + sp[2] / 2, step=sp[2], length=sz[2])
    xs, ys
  elseif nd == 3
    xs = range(or[1] + sp[1] / 2, step=sp[1], length=sz[1])
    ys = range(or[2] + sp[2] / 2, step=sp[2], length=sz[2])
    zs = range(or[3] + sp[3] / 2, step=sp[3], length=sz[3])
    xs, ys, zs
  else
    throw(ErrorException("not implemented"))
  end
end
