# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{TransformedGrid}})
  # add colorant to compute graph
  colorant!(plot)

  # add parent domain and transformation to compute graph
  Makie.map!(plot, :object, [:pdom, :trans]) do tdom
    parent(tdom), Meshes.transform(tdom)
  end

  if isoptimized(plot.trans[])
    # visualize parent domain and transform visualization
    viz!(plot, plot.attributes, plot.pdom)
    makietransform!(plot, plot.trans[])
  else
    # visualize as grid
    vizgrid!(plot)
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
