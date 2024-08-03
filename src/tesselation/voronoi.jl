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
  rawval = map(p -> CoordRefSystems.raw(coords(p)), pset)
  triang = triangulate(rawval, rng=method.rng)
  vorono = voronoi(triang, clip=true)

  # mesh with all (possibly unused) points
  points = map(get_polygon_points(vorono)) do xy
    coords = CoordRefSystems.reconstruct(C, T.(xy))
    Point(coords)
  end
  polygs = get_polygons(vorono)
  connec = Vector{Connectivity}(undef, length(polygs))
  for (order_idx, vertices_inds) in polygs
    tup = ntuple(i -> vertices_inds[i], length(vertices_inds) - 1)
    connec[order_idx] = connect(tup, Ngon)
  end
  mesh = SimpleMesh(points, connec)

  # remove unused points
  mesh |> Repair{1}()
end
