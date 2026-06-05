# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function vizgrid!(plot::Viz{<:Tuple{TransformedGrid}}, M::Type{<:𝔼}, pdim::Val, edim::Val)
  # parent grid and transformation
  Makie.map!(plot, :object, [:pgrid, :trans]) do grid
    pgrid = parent(grid)
    trans = transformation(grid)
    pgrid, trans
  end

  if isoptimized(plot.trans[])
    # visualize parent grid and transform visualization
    plot.object = plot.pgrid
    vizgrid!(plot)
    makietransform!(plot)
  else
    # fallback to full grid visualization
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

makietransform!(plot) = makietransform!(plot, plot.trans[])

function makietransform!(plot, trans::Rotate{<:Angle2d})
  rot = first(TB.parameters(trans))
  θ = first(Rotations.params(rot))
  Makie.rotate!(plot, θ)
end

function makietransform!(plot, trans::Translate)
  offsets = first(TB.parameters(trans))
  Makie.translate!(plot, ustrip.(offsets)...)
end

function makietransform!(plot, trans::Scale)
  factors = first(TB.parameters(trans))
  Makie.scale!(plot, factors...)
end

function makietransform!(plot, trans::Affine{2})
  A, b = TB.parameters(trans)
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

makietransform!(plot, trans::TB.Identity) = nothing

makietransform!(plot, trans::TB.SequentialTransform) = foreach(t -> makietransform!(plot, t), trans)
