# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

# preprocess colors provided by user
process(values::AbstractVector, colorscheme, colorrange, alphas) = colorfy(values; alphas, colorscheme, colorrange)
process(value, colorscheme, colorrange, alphas) = process([value], colorscheme, colorrange, alphas) |> first

# assumes that meshes and colors have the same length
function vizmany!(plot, meshes, colors)
  meshid = Symbol(meshes, :_meshes)
  colorid = Symbol(meshes, :_color)
  Makie.map!(plot, [meshes, colors], [meshid, colorid]) do ms, cs
    mesh = reduce(merge, ms)
    color = [cs[i] for (i, m) in enumerate(ms) for _ in 1:nelements(m)]
    mesh, color
  end
  viz!(plot, plot[meshid], color=plot[colorid], pointsize=plot.pointsize, segmentsize=plot.segmentsize)
end

asray(v::Vec) = Ray(Point(zero(v)...), v)

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
