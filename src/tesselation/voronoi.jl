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

### Notes

Wraps DelaunayTriangulation.jl. For any internal errors, file an issue at 
[DelaunayTriangulation.jl](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/issues/new)
"""
struct VoronoiTesselation{RNG<:AbstractRNG} <: TesselationMethod
  rng::RNG
end

VoronoiTesselation(rng=Random.default_rng()) = VoronoiTesselation(rng)

function tesselate(pset::PointSet, method::VoronoiTesselation)
  C = crs(pset)
  T = numtype(lentype(pset))
  assertion(CoordRefSystems.ncoords(C) == 2, "points must have 2 coordinates")

  # perform tesselation with raw coordinates
  rawval = map(p -> CoordRefSystems.rawvalues(coords(p)), pset)
  triang = triangulate(rawval, rng=method.rng)
  vorono = voronoi(triang, clip=true)

  # mesh with all (possibly unused) points
  points = map(get_polygon_points(vorono)) do xy
    coords = CoordRefSystems.reconstruct(C, T.(xy))
    Point(coords)
  end
  polygs = each_polygon(vorono)
  tuples = [Tuple(inds[begin:(end - 1)]) for inds in polygs]
  connec = connect.(tuples)
  mesh = SimpleMesh(points, connec)

  # remove unused points
  mesh |> Repair{1}()
end
