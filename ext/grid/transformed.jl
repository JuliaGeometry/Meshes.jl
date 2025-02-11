# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function vizgrid!(plot::Viz{<:Tuple{TransformedGrid}}, M::Type{<:𝔼}, pdim::Val, edim::Val)
  tgrid = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colormap = plot[:colormap]
  colorrange = plot[:colorrange]
  showsegments = plot[:showsegments]
  segmentcolor = plot[:segmentcolor]
  segmentsize = plot[:segmentsize]

  # retrieve transformation
  trans = Makie.@lift Meshes.transform($tgrid)

  if isoptimized(trans[]) # visualize parent grid and transform visualization
    grid = Makie.@lift parent($tgrid)
    viz!(plot, grid; color, alpha, colormap, showsegments, segmentcolor, segmentsize)
    makietransform!(plot, trans)
  elseif pdim == Val(2) # visualize quadrangle mesh with texture using uv coords
    # decide whether or not to reverse connectivity list
    rfunc = Makie.@lift _reverse(crs($tgrid))

    verts = Makie.@lift map(asmakie, vertices($tgrid))
    quads = Makie.@lift [GB.QuadFace($rfunc(indices(e))) for e in elements(topology($tgrid))]

    colorant = Makie.@lift process($color, $colormap, $colorrange, $alpha)

    nverts = Makie.@lift length($verts)
    nquads = Makie.@lift length($quads)
    ncolor = Makie.@lift length($colorant)

    dims = Makie.@lift size($tgrid)
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

# --------------
# OPTIMIZATIONS
# --------------

isoptimized(t) = false
isoptimized(::Rotate{<:Angle2d}) = true
isoptimized(::Translate) = true
isoptimized(::Scale) = true
function isoptimized(t::Affine{2})
  A, _ = TB.parameters(t)
  isdiag(A) || isrotation(A)
end
isoptimized(::TB.Identity) = true
isoptimized(t::TB.SequentialTransform) = all(isoptimized, t)

function makietransform!(plot, trans::Makie.Observable{<:Rotate{<:Angle2d}})
  rot = first(TB.parameters(trans[]))
  θ = first(Rotations.params(rot))
  Makie.rotate!(plot, θ)
end

function makietransform!(plot, trans::Makie.Observable{<:Translate})
  offsets = first(TB.parameters(trans[]))
  Makie.translate!(plot, ustrip.(offsets)...)
end

function makietransform!(plot, trans::Makie.Observable{<:Scale})
  factors = first(TB.parameters(trans[]))
  Makie.scale!(plot, factors...)
end

function makietransform!(plot, trans::Makie.Observable{<:Affine{2}})
  A, b = TB.parameters(trans[])
  if isdiag(A)
    s₁, s₂ = A[1, 1], A[2, 2]
    Makie.scale!(plot, s₁, s₂)
  else
    rot = convert(Angle2d, A)
    θ = first(Rotations.params(rot))
    Makie.rotate!(plot, θ)
  end
  Makie.translate!(plot, ustrip.(b)...)
end

makietransform!(plot, trans::Makie.Observable{<:TB.Identity}) = nothing

makietransform!(plot, trans::Makie.Observable{<:TB.SequentialTransform}) =
  foreach(t -> makietransform!(plot, Makie.Observable(t)), trans[])
