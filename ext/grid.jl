# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

Makie.plot!(plot::Viz{<:Tuple{Grid}}) = vizgrid!(plot)

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
  showsegments = plot.showsegments

  # process color spec into colorant
  Makie.map!(process, plot, [:color, :colormap, :colorrange, :alpha], :colorant)

  # number of vertices, elements and colors
  Makie.map!(nvertices, plot, [:object], :nverts)
  Makie.map!(nelements, plot, [:object], :nelems)
  Makie.map!(plot, [:colorant], :ncolor) do colorant
    colorant isa AbstractVector ? length(colorant) : 1
  end

  # visualize quadrangle mesh with texture using uv coords
  # plots with uv coords are always interpolated,
  # so it is only used in the case ncolor == nverts
  # or when there is a large number of elements
  if pdim === Val(2) && (plot.ncolor[] == 1 || plot.ncolor[] == plot.nverts[] || plot.nelems[] ≥ 1000)
    # decide whether or not to reverse connectivity list
    Makie.map!(plot, [:object], :rfunc) do grid
      crs(grid) <: LatLon && orientation(first(grid)) == CW ? reverse : identity
    end

    Makie.map!(plot, [:object], :verts) do grid
      map(asmakie, eachvertex(grid))
    end
    Makie.map!(plot, [:object, :rfunc], :quads) do grid, rfunc
      [GB.QuadFace(rfunc(indices(e))) for e in elements(topology(grid))]
    end

    Makie.map!(size, plot, [:object], :dims)
    Makie.map!(Meshes.vsize, plot, [:object], :vdims)

    Makie.map!(
      plot,
      [:colorant, :ncolor, :nelems, :nverts, :dims, :vdims],
      :texture
    ) do colorant, ncolor, nelems, nverts, dims, vdims
      if ncolor == 1
        fill(colorant, dims)
      elseif ncolor == nelems
        reshape(colorant, dims)
      elseif ncolor == nverts
        reshape(colorant, vdims)
      else
        throw(ArgumentError("invalid number of colors"))
      end
    end

    Makie.map!(plot, [:vdims], :uv) do vdims
      [Makie.Vec2f(v, 1 - u) for v in range(0, 1, vdims[2]) for u in range(0, 1, vdims[1])]
    end

    Makie.map!(plot, [:verts, :quads, :uv], :mesh) do verts, quads, uv
      GB.Mesh(verts, quads, uv=uv)
    end

    # enable shading in 3D
    shading = edim === Val(3)

    Makie.mesh!(plot, plot.mesh, color=plot.texture, shading=shading)

    if showsegments[]
      vizfacets!(plot)
    end
  else # fallback to triangle mesh visualization
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
include("grid/transformed.jl")

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
