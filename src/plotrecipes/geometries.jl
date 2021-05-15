# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

@recipe function f(ray::Ray)
  seriestype --> :path
  linecolor --> :black
  arrow --> true
  [ray(0), ray(1)]
end

@recipe function f(segment::Segment)
  seriestype --> :path
  linecolor --> :black
  vertices(segment)
end

@recipe function f(sphere::Sphere, nsamples=100)
  seriestype --> :path
  linecolor --> :black
  samples = sample(sphere, RegularSampling(nsamples))
  points  = collect(samples)
  [points; first(points)]
end

@recipe function f(ball::Ball, nsamples=100)
  seriestype --> :shape
  linecolor --> :black
  samples = sample(boundary(ball), RegularSampling(nsamples))
  points  = collect(samples)
  [points; first(points)]
end

@recipe function f(polygon::Polygon)
  seriestype --> :path
  seriescolor --> :auto
  fill --> true
  points = vertices(polygon)
  [points; first(points)]
end

@recipe function f(chain::Chain)
  seriestype --> :path
  linecolor --> :black
  points = vertices(chain)
  isclosed(chain) ? [points; first(points)] : points
end

@recipe function f(polyarea::PolyArea)
  if hasholes(polyarea)
    mesh  = discretize(polyarea, FIST())
    geoms = collect(mesh)
  else
    geoms = [first(chains(polyarea))]
  end

  seriestype --> :path
  @series begin
    geoms
  end
  for chain in chains(polyarea)
    @series begin
      primary --> false
      chain
    end
  end
end

@recipe function f(multi::Multi)
  collect(multi)
end

@recipe function f(geometries::AbstractVector{G}) where {G<:Geometry}
  @series begin
    first(geometries)
  end
  for geometry in Iterators.drop(geometries, 1)
    @series begin
      primary --> false
      geometry
    end
  end
end
