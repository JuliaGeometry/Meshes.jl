# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{Grid}})
  grid = plot[:object]
  M = Makie.@lift manifold($grid)
  pdim = Makie.@lift paramdim($grid)
  edim = Makie.@lift embeddim($grid)
  vizgrid!(plot, M[], Val(pdim[]), Val(edim[]))
end

function vizgrid!(plot, ::Type{<:🌐}, pdim::Val, edim::Val)
  vizgrid!(plot, 𝔼, pdim, edim)
end

vizgrid!(plot, M::Type{<:𝔼}, pdim::Val, edim::Val) = vizgridfallback!(plot, M, pdim, edim)

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

function vizgridfallback!(plot, M, pdim, edim)
  grid = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colormap = plot[:colormap]
  colorrange = plot[:colorrange]

  if pdim == Val(2) # visualize quadrangle mesh with texture using uv coords
    # decide whether or not to reverse connectivity list
    rfunc = Makie.@lift _reverse(crs($grid))

    verts = Makie.@lift map(asmakie, vertices($grid))
    quads = Makie.@lift [GB.QuadFace($rfunc(indices(e))) for e in elements(topology($grid))]

    colorant = Makie.@lift process($color, $colormap, $colorrange, $alpha)

    nverts = Makie.@lift length($verts)
    nquads = Makie.@lift length($quads)
    ncolor = Makie.@lift $colorant isa AbstractVector ? length($colorant) : 1

    dims = Makie.@lift size($grid)
    texture = if ncolor[] == 1
      Makie.@lift fill($colorant, $dims)
    elseif ncolor[] == nquads[]
      Makie.@lift reshape($colorant, $dims)
    elseif ncolor[] == nverts[]
      Makie.@lift reshape($colorant, $dims .+ 1)
    else
      throw(ArgumentError("invalid number of colors"))
    end

    uv = Makie.@lift [Makie.Vec2f(v, 1 - u) for v in range(0, 1, $dims[2] + 1) for u in range(0, 1, $dims[1] + 1)]

    mesh = Makie.@lift GB.Mesh(Makie.meta($verts, uv=$uv), $quads)

    shading = edim == Val(3) ? Makie.FastShading : Makie.NoShading

    Makie.mesh!(plot, mesh, color=texture, shading=shading)
  else # fallback to triangle mesh visualization
    vizmesh!(plot, M, pdim, edim)
  end
end

_reverse(::Type{<:CRS}) = identity
_reverse(::Type{<:LatLon}) = reverse

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
