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

function tesselate(pset::PointSet, method::VoronoiTesselation)
  C = crs(pset)
  assertion(CoordRefSystems.ncoords(C) == 2, "the number of coordinates of the points must be 2")

  # perform tesselation with raw coordinates
  cs = map(p -> CoordRefSystems.rawvalues(coords(p)), pset)
  triang = triangulate(cs, rng=method.rng)
  vorono = voronoi(triang, clip=true)

  # mesh with all (possibly unused) points
  points = map(cs -> Point(CoordRefSystems.reconstruct(C, cs)), get_polygon_points(vorono))
  polygs = each_polygon(vorono)
  tuples = [Tuple(inds[begin:(end - 1)]) for inds in polygs]
  connec = connect.(tuples)
  mesh = SimpleMesh(points, connec)

  # remove unused points
  mesh |> Repair{1}()
end
