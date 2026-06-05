# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{SubDomain}})
  # retrieve parent domain
  Makie.map!(parent, plot, :object, :pdom)

  if plot.pdom[] isa CartesianGrid
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
  Makie.map!(plot, :object, [:coords, :marker, :shading]) do sgrid
    pgrid = parent(sgrid)
    nd = embeddim(pgrid)
    sp = ustrip.(spacing(pgrid))

    # coordinates of markers
    coords = map(sgrid) do e
      ustrip.(to(centroid(e))) .+ sp ./ 2
    end

    # rectangle markers
    marker = Makie.Rect{length(sp)}(-1 .* sp, sp)

    # enable shading in 3D
    shading = nd == 3

    coords, marker, shading
  end

  # all geometries are equal, use mesh scatter
  Makie.meshscatter!(plot, plot.coords, marker=plot.marker, markersize=1, color=plot.colorant, shading=plot.shading)
end
