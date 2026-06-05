# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function Makie.plot!(plot::Viz{<:Tuple{SubDomain}})
  Makie.map!(process, plot, [:color, :colormap, :colorrange, :alpha], :colorant)
  vizsubdom!(plot)
end

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

vizsubdom!(plot, M::Type, pdim::Val, edim::Val) = vizsubdomfallback!(plot, M, pdim, edim)

function vizsubdomfallback!(plot, ::Type, ::Val, ::Val)
  # visualize as geometry set
  gset = GeometrySet(collect(plot.object[]))
  Mke.update!(plot, object=gset)
  vizgset!(plot)
end

# ----------------
# SPECIALIZATIONS
# ----------------

const SubCartesianGrid{M,CRS} = SubDomain{M,CRS,<:CartesianGrid}

function vizsubdom!(plot::Viz{<:Tuple{SubCartesianGrid}}, ::Type, ::Val, edim::Val)
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
