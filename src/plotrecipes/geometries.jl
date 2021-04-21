# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@recipe function f(geometries::AbstractVector{<:Geometry})
  for geometry in geometries
    @series begin
      geometry
    end
  end
end

@recipe function f(multi::Multi)
  @series begin
    first(multi)
  end
  for geometry in Iterators.drop(multi, 1)
    @series begin
      primary --> false
      geometry
    end
  end
end

@recipe function f(segment::Segment)
  seriestype --> :path
  linecolor --> :black

  vertices(segment)
end

@recipe function f(polygon::Polygon)
  seriestype --> :path
  seriescolor --> :auto

  points = vertices(polygon)
  [points; first(points)]
end

@recipe function f(ray::Ray)
  seriestype --> :path
  linecolor --> :black
  arrow --> true

  [ray(0), ray(1)]
end

@recipe function f(sphere::Sphere, nsamples=100)
  seriestype --> :path
  linecolor --> :black

  samples = sample(sphere, RegularSampling(nsamples))
  points  = collect(samples)
  [points; first(points)]
end

@recipe function f(chain::Chain)
  seriestype --> :path
  linecolor --> :black

  points = vertices(chain)

  isclosed(chain) ? [points; first(points)] : points
end

@recipe function f(polyarea::PolyArea)
  seriestype --> :path
  seriescolor --> :auto
  linecolor --> :black
  fill --> true

  if hasholes(polyarea)
    mesh  = discretize(polyarea, FIST())
    geoms = collect(mesh)
  else
    geoms = [first(chains(polyarea))]
  end

  @series begin
    first(geoms)
  end

  for geom in Iterators.drop(geoms, 1)
    @series begin
      primary --> false
      geom
    end
  end
end
