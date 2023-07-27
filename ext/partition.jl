# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

Makie.plottype(::Partition) = Viz{<:Tuple{Partition}}

function Makie.plot!(plot::Viz{<:Tuple{Partition}})
  # retrieve partition object
  partition = plot[:object][]

  colors = distinguishable_colors(length(partition), transform=protanopic)

  for (i, subset) in Iterators.enumerate(partition)
    # fallback to collection recipe
    viz!(plot, subset,
      color       = colors[i],
      alpha       = plot[:alpha],
      facetcolor  = plot[:facetcolor],
      showfacets  = plot[:showfacets],
      pointsize   = plot[:pointsize],
      segmentsize = plot[:segmentsize]
    )
  end
end
