# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# assumes that meshes and colors have the same length
# and that colors are processed with alpha and colormap
function vizmany!(plot, meshes, colors)
  pointsize = plot.pointsize
  segmentsize = plot.segmentsize

  rmesh = Symbol(meshes, :_rmeshes)
  rcolor = Symbol(meshes, :_rcolor)

  Makie.map!(plot, [meshes, colors], [rmesh, rcolor]) do ms, cs
    rmes = reduce(merge, ms)
    color = [cs[i] for (i, m) in enumerate(ms) for _ in 1:nelements(m)]
    (rmes, color)
  end

  viz!(plot, plot[rmesh], color=plot[rcolor], pointsize=pointsize, segmentsize=segmentsize)
end

asray(vec::Vec{Dim,ℒ}) where {Dim,ℒ} = Ray(Point(ntuple(i -> zero(ℒ), Dim)), vec)

asmakie(v::Vec) = Makie.Vec{length(v),numtype(eltype(v))}(Tuple(ustrip.(v)))

asmakie(p::Point) = Makie.Point{embeddim(p),numtype(lentype(p))}(Tuple(ustrip.(to(p))))

function asmakie(poly::Polygon)
  rs = rings(poly)
  outer = map(asmakie, eachvertex(rs[1]))
  if hasholes(poly)
    inners = map(i -> map(asmakie, eachvertex(rs[i])), 2:length(rs))
    Makie.Polygon(outer, inners)
  else
    Makie.Polygon(outer)
  end
end
