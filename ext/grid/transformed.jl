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
    viz!(plot, grid; color, alpha, colormap, colorrange, showsegments, segmentcolor, segmentsize)
    makietransform!(plot, trans)
  else # fallback to full grid visualization
    vizgridfallback!(plot, M, pdim, edim)
  end
end

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
