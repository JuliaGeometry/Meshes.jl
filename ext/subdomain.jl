# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{SubDomain}})
  if parent(plot.object[]) isa CartesianGrid
    # visualize as subset of Cartesian grid
    colorant!(plot)
    vizsubcartgrid!(plot)
  else
    # visualize as geometry set
    Makie.map!(sdom -> GeometrySet(collect(sdom)), plot, :object, :gset)
    viz!(
      plot,
      plot.gset,
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
  end
end

# --------------
# OPTIMIZATIONS
# --------------

function vizsubcartgrid!(plot)
  # retrieve grid paramaters
  Makie.map!(plot, :object, [:xyz, :rec]) do sgrid
    sp = ustrip.(spacing(parent(sgrid)))

    # coordinates of markers
    xyz = map(sgrid) do e
      ustrip.(to(centroid(e))) .+ sp ./ 2
    end

    # rectangle markers
    rec = Makie.Rect{length(sp)}(-1 .* sp, sp)

    xyz, rec
  end

  # enable shading in 3D
  shading = embeddim(plot.object[]) == 3

  # all geometries are equal, use mesh scatter
  Makie.meshscatter!(plot, plot.xyz, marker=plot.rec, markersize=1, color=plot.colorant, shading=shading)
end
