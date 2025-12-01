# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

function vizmany!(plot, meshes, colors)
  alpha = plot[:alpha]
  colormap = plot[:colormap]
  pointsize = plot[:pointsize]
  segmentsize = plot[:segmentsize]

  mesh = Makie.@lift reduce(merge, $meshes)
  color = Makie.@lift [$colors[i] for (i, m) in enumerate($meshes) for _ in 1:nelements(m)]

  viz!(plot, mesh, color=color, alpha=alpha, colormap=colormap, pointsize=pointsize, segmentsize=segmentsize)
end

asray(vec::Vec{Dim,ℒ}) where {Dim,ℒ} = Ray(Point(ntuple(i -> zero(ℒ), Dim)), vec)

asmakie(v::Vec) = Makie.Vec{length(v),numtype(eltype(v))}(ustrip.(Tuple(v)))

asmakie(p::Point) = Makie.Point{embeddim(p),numtype(lentype(p))}(ustrip.(Tuple(to(p))))

asmakie(b::Box) = Makie.Rect([asmakie(p) for p in boundarypoints(b)])

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
