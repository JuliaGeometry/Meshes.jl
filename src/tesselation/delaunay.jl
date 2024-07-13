# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DelaunayTesselation([rng])

Unconstrained Delaunay tesselation of point sets.
Optionally, specify the random number generator `rng`.

## References

* Cheng et al. 2012. [Delaunay Mesh Generation]
  (https://people.eecs.berkeley.edu/~jrs/meshbook.html)

### Notes

Wraps DelaunayTriangulation.jl. For any internal errors, file an issue at 
[DelaunayTriangulation.jl](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/issues/new)
"""
struct DelaunayTesselation{RNG<:AbstractRNG} <: TesselationMethod
  rng::RNG
end

DelaunayTesselation(rng=Random.default_rng()) = DelaunayTesselation(rng)

function tesselate(pset::PointSet, method::DelaunayTesselation)
  assertion(CoordRefSystems.ncoords(crs(pset)) == 2, "points must have 2 coordinates")

  # perform tesselation with raw coordinates
  rawval = map(p -> CoordRefSystems.rawvalues(coords(p)), pset)
  triang = triangulate(rawval, rng=method.rng)
  connec = connect.(each_solid_triangle(triang))
  SimpleMesh(collect(pset), connec)
end
