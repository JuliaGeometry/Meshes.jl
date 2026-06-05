# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{TransformedDomain}})
  # retrieve parent domain and transformation
  Makie.map!(plot, :object, [:pdom, :trans]) do tdom
    pdom = parent(tdom)
    trans = transformation(tdom)
    pdom, trans
  end

  if isoptimized(plot.trans[])
    # visualize parent domain and transform visualization
    viz!(
      plot,
      plot.pdom,
      color=plot.color,
      alpha=plot.alpha,
      colormap=plot.colormap,
      colorrange=plot.colorrange,
      showsegments=plot.showsegments,
      segmentcolor=plot.segmentcolor,
      segmentsize=plot.segmentsize,
      showpoints=plot.showpoints,
      pointmarker=plot.pointmarker,
      pointcolor=plot.pointcolor,
      pointsize=plot.pointsize
    )
    makietransform!(plot, plot.trans[])
  elseif plot.pdom[] isa Grid
    # visualize as grid
    colorant!(plot)
    vizgrid!(plot)
  else
    error("not implemented")
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
