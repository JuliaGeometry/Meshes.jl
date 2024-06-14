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
"""
struct DelaunayTesselation{RNG<:AbstractRNG} <: TesselationMethod
  rng::RNG
end

DelaunayTesselation(rng=Random.default_rng()) = DelaunayTesselation(rng)

function tesselate(pset::PointSet{2}, method::DelaunayTesselation)
  # perform tesselation with raw coordinates
  coords = map(p -> ustrip.(to(p)), pset)
  triang = triangulate(coords, rng=method.rng)
  connec = connect.(each_solid_triangle(triang))
  SimpleMesh(collect(pset), connec)
end
