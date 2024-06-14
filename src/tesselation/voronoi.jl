# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    VoronoiTesselation([rng])

Unconstrained Voronoi tesselation of point sets.
Optionally, specify the random number generator `rng`.

## References

* Cheng et al. 2012. [Delaunay Mesh Generation]
  (https://people.eecs.berkeley.edu/~jrs/meshbook.html)
"""
struct VoronoiTesselation{RNG<:AbstractRNG} <: TesselationMethod
  rng::RNG
end

VoronoiTesselation(rng=Random.default_rng()) = VoronoiTesselation(rng)

function tesselate(pset::PointSet{2}, method::VoronoiTesselation)
  # perform tesselation with raw coordinates
  coords = map(p -> ustrip.(to(p)), pset)
  triang = triangulate(coords, rng=method.rng)
  vorono = voronoi(triang, clip=true)

  # mesh with all (possibly unused) points
  points = get_polygon_points(vorono)
  polygs = each_polygon(vorono)
  tuples = [Tuple(inds[begin:(end - 1)]) for inds in polygs]
  connec = connect.(tuples)
  mesh = SimpleMesh(points, connec)

  # remove unused points
  mesh |> Repair{1}()
end
