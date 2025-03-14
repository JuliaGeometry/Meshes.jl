# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

mayberepeat(value, objs) = value

mayberepeat(value::AbstractVector, objs) = [value[e] for (e, obj) in enumerate(objs) for _ in 1:length(obj)]

mayberepeat(value::AbstractVector, objs::AbstractVector{<:Multi}) =
  [value[e] for (e, obj) in enumerate(objs) for _ in 1:length(parent(obj))]

mayberepeat(value::AbstractVector, objs::AbstractVector{<:Geometry}) = value

concat(obj₁, obj₂) = vcat(obj₁, obj₂)

concat(mesh₁::Mesh, mesh₂::Mesh) = merge(mesh₁, mesh₂)

function vizmany!(plot, objs, color)
  alpha = plot[:alpha]
  colormap = plot[:colormap]
  pointsize = plot[:pointsize]
  segmentsize = plot[:segmentsize]

  object = Makie.@lift reduce(concat, $objs)
  colors = Makie.@lift mayberepeat($color, $objs)
  alphas = Makie.@lift mayberepeat($alpha, $objs)

  viz!(plot, object, color=colors, alpha=alphas, colormap=colormap, pointsize=pointsize, segmentsize=segmentsize)
end

asmakie(geoms::AbstractVector{<:Geometry}) = asmakie.(geoms)

asmakie(multis::AbstractVector{<:Multi}) = mapreduce(m -> asmakie.(parent(m)), vcat, multis)

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

asmakie(p::Point) = Makie.Point{embeddim(p),numtype(lentype(p))}(ustrip.(Tuple(to(p))))

asmakie(v::Vec) = Makie.Vec{length(v),numtype(eltype(v))}(ustrip.(Tuple(v)))

_discretize(geom) = discretize(geom)

function _discretize(box::Box{🌐})
  T = numtype(Meshes.lentype(box))
  discretize(box, MaxLengthDiscretization(T(100) * u"km"))
end

function _discretize(chain::Chain{🌐})
  T = numtype(Meshes.lentype(chain))
  discretize(chain, MaxLengthDiscretization(T(1000) * u"km"))
end
