# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{SubDomain}})
  # collect geometries into geometry set
  Makie.map!(plot, :object, :gset) do subdom
    GeometrySet(collect(subdom))
  end

  # visualize as geometry set
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

# ----------------
# SPECIALIZATIONS
# ----------------

const SubCartesianGrid{M,CRS} = SubDomain{M,CRS,<:CartesianGrid}

function Makie.plot!(plot::Viz{<:Tuple{SubCartesianGrid}})
  colorant!(plot)

  # retrieve grid paramaters
  Makie.map!(plot, :object, [:scoords, :smarker]) do subgrid
    grid = parent(subgrid)
    sp = ustrip.(spacing(grid))

    # coordinates of markers
    scoords = map(subgrid) do e
      ustrip.(to(centroid(e))) .+ sp ./ 2
    end

    # rectangle markers
    smarker = Makie.Rect{length(sp)}(-1 .* sp, sp)

    scoords, smarker
  end

  # enable shading in 3D
  shading = edim == Val(3)

  # all geometries are equal, use mesh scatter
  Makie.meshscatter!(plot, plot.scoords, marker=plot.smarker, markersize=1, color=plot.colorant, shading=shading)
end
