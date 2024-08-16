# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function vizgrid!(plot::Viz{<:Tuple{Meshes.TransformedGrid}}, M::Type{<:𝔼}, pdim::Val, edim::Val)
  tgrid = plot[:object]
  color = plot[:color]
  alpha = plot[:alpha]
  colormap = plot[:colormap]
  colorrange = plot[:colorrange]
  grid = Makie.@lift parent($tgrid)
  trans = Makie.@lift Meshes.transform($tgrid)
  if isoptimized(crs(grid[]), trans[])
    showsegments = plot[:showsegments]
    segmentcolor = plot[:segmentcolor]
    segmentsize = plot[:segmentsize]
    viz!(plot, grid; color, alpha, colormap, showsegments, segmentcolor, segmentsize)
    makietransform!(plot, trans)
  else
    if paramdim(tgrid[]) == 2
      colorant = Makie.@lift process($color, $colormap, $colorrange, $alpha)
      texture = Makie.@lift reshape($colorant, size($tgrid))
      coords = Makie.@lift map(asmakie, vertices($tgrid))
      quads = Makie.@lift [GB.QuadFace(indices(e)) for e in elements(topology($tgrid))]
      sz = Makie.@lift size($tgrid)
      uv = Makie.@lift [Makie.Vec2f(u, v) for u in range(0, 1, $sz[1]) for v in range(0, 1, $sz[2])]
      msh = Makie.@lift GB.Mesh(Makie.meta($coords, uv=$uv), $quads)
      Makie.mesh!(msh, color=texture)
    else
      vizmesh!(plot, M, pdim, edim)
    end
  end
end

isoptimized(::Type, ::TB.Identity) = true
isoptimized(CRS::Type, t::TB.SequentialTransform) = all(tᵢ -> isoptimized(CRS, tᵢ), t)

isoptimized(::Type, ::GeometricTransform) = false

isoptimized(::Type{<:Cartesian2D}, ::Proj{<:Projected}) = true
isoptimized(::Type{<:Projected}, ::Proj{<:Cartesian2D}) = true

isoptimized(::Type, ::Rotate{<:Angle2d}) = true
isoptimized(::Type, ::Translate) = true
isoptimized(::Type, ::Scale) = true
function isoptimized(::Type, t::Affine{2})
  A, _ = TB.parameters(t)
  isdiag(A) || isrotation(A)
end

makietransform!(plot, trans::Makie.Observable{<:TB.Identity}) = nothing

makietransform!(plot, trans::Makie.Observable{<:Proj}) = nothing

makietransform!(plot, trans::Makie.Observable{<:TB.SequentialTransform}) =
  foreach(t -> makietransform!(plot, Makie.Observable(t)), trans[])

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
