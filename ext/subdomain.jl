# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

Makie.plot!(plot::Viz{<:Tuple{SubDomain}}) = vizsubdom!(plot)

function vizsubdom!(plot)
  subdom = plot.object[]
  M = manifold(subdom)
  pdim = paramdim(subdom)
  edim = embeddim(subdom)
  vizsubdom!(plot, M, Val(pdim), Val(edim))
end

# ---------------
# IMPLEMENTATION
# ---------------

function vizsubdom!(plot, ::Type{<:🌐}, pdim::Val, edim::Val)
  vizsubdom!(plot, 𝔼, pdim, edim)
end

function vizsubdom!(plot, ::Type{<:𝔼}, ::Val, ::Val)
  # construct geometry set
  Makie.map!(plot, :object, :gset) do subdom
    GeometrySet(collect(subdom))
  end

  # forward attributes
  viz!(
    plot,
    plot.gset;
    plot.color,
    plot.alpha,
    plot.colormap,
    plot.colorrange,
    plot.showsegments,
    plot.segmentcolor,
    plot.segmentsize,
    plot.showpoints,
    plot.pointmarker,
    plot.pointcolor,
    plot.pointsize
  )
end

const SubCartesianGrid{M,CRS} = SubDomain{M,CRS,<:CartesianGrid}

function vizsubdom!(plot::Viz{<:Tuple{SubCartesianGrid}}, ::Type{<:𝔼}, ::Val, edim::Val)
  # process color spec into colorant
  Makie.map!(process, plot, [:color, :colormap, :colorrange, :alpha], :colorant)

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
