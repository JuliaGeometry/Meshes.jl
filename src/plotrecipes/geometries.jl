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

@recipe function f(segment::Segment)
  seriestype --> :path
  seriescolor --> :black
  primary --> false

  vertices(segment)
end

@recipe function f(polygon::Polygon)
  seriestype --> :path
  seriescolor --> :black
  primary --> false

  points = vertices(polygon)
  [points; first(points)]
end

@recipe function f(ray::Ray)
  seriestype --> :path
  seriescolor --> :black
  arrow --> true
  label --> "ray"

  [ray(0), ray(1)]
end

@recipe function f(sphere::Sphere, nsamples=100)
  seriestype --> :path
  seriescolor --> :black
  label --> "sphere"

  samples = sample(sphere, RegularSampling(nsamples))
  points  = collect(samples)
  [points; first(points)]
end

@recipe function f(chain::Chain)
  seriestype --> :path
  seriescolor --> :black
  label --> "chain"

  points = vertices(chain)

  isclosed(chain) ? [points; first(points)] : points
end

@recipe function f(polyarea::PolyArea)
  seriestype --> :path
  seriescolor --> :black
  label --> "polyarea"

  pchains = chains(polyarea)

  # plot outer chain
  @series begin
    first(pchains)
  end

  # plot inner chains
  for pchain in pchains[2:end]
    @series begin
      primary --> false
      pchain
    end
  end
end
