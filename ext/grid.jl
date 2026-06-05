# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{Grid}})
  Makie.map!(process, plot, [:color, :colormap, :colorrange, :alpha], :colorant)
  vizgrid!(plot)
end

function vizgrid!(plot)
  grid = plot.object[]
  M = manifold(grid)
  pdim = paramdim(grid)
  edim = embeddim(grid)
  vizgrid!(plot, M, Val(pdim), Val(edim))
end

# ---------------
# IMPLEMENTATION
# ---------------

vizgrid!(plot, M::Type, pdim::Val, edim::Val) = vizgridfallback!(plot, M, pdim, edim)

function vizgridfallback!(plot, M, pdim, edim)
  # number of vertices, elements and colors
  Makie.map!(plot, [:object, :colorant], [:nv, :ne, :nc]) do grid, colorant
    nv = nvertices(grid)
    ne = nelements(grid)
    nc = colorant isa AbstractVector ? length(colorant) : 1
    nv, ne, nc
  end

  # plots with uv coords are always interpolated, so it is
  # only used in the case nc == nv or when the grid is large
  if pdim === Val(2) && (plot.nc[] == 1 || plot.nc[] == plot.nv[] || plot.ne[] ≥ 1000)
    # visualize quadrangle mesh with texture using uv coords
    Makie.map!(plot, [:object, :colorant], [:mesh, :texture]) do grid, colorant
      # retrieve relevant sizes
      sz = size(grid)
      vz = Meshes.vsize(grid)

      # decide whether or not to reverse connectivity list
      rev = crs(grid) <: LatLon && orientation(first(grid)) == CW ? reverse : identity

      # vertices and quadrangles
      verts = map(asmakie, eachvertex(grid))
      quads = [GB.QuadFace(rev(indices(e))) for e in elements(topology(grid))]

      # uv coordinates for texture mapping
      uv = [Makie.Vec2f(v, 1 - u) for v in range(0, 1, vz[2]) for u in range(0, 1, vz[1])]

      mesh = GB.Mesh(verts, quads; uv)

      # texture map
      texture = if plot.nc[] == 1
        fill(colorant, sz)
      elseif plot.nc[] == plot.nv[]
        reshape(colorant, vz)
      elseif plot.nc[] == plot.ne[]
        reshape(colorant, sz)
      else
        throw(ArgumentError("invalid number of colors"))
      end

      mesh, texture
    end

    # enable shading in 3D
    shading = edim === Val(3)

    Makie.mesh!(plot, plot.mesh, color=plot.texture, shading=shading)

    if plot.showsegments[]
      vizfacets!(plot)
    end
  else
    # fallback to triangle mesh visualization
    vizmesh!(plot)
  end
end

# -------
# FACETS
# -------

function vizfacets!(plot::Viz{<:Tuple{Grid}})
  grid = plot.object[]
  M = manifold(grid)
  pdim = paramdim(grid)
  edim = embeddim(grid)
  vizgridfacets!(plot, M, Val(pdim), Val(edim))
end

vizgridfacets!(plot, M::Type, pdim::Val, edim::Val) = vizmeshfacets!(plot, M, pdim, edim)

# ----------------
# SPECIALIZATIONS
# ----------------

include("grid/cartesian.jl")
include("grid/rectilinear.jl")
include("grid/structured.jl")

# -----------------
# HELPER FUNCTIONS
# -----------------

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
