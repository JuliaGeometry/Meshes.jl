# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DelaunayTriangulation()

Constrained Delaunay triangulation of polygons.
Optionally, specify the random number generator `rng`.

!!! note "Wrapper"

    Wraps DelaunayTriangulation.jl's `triangulate`. For any internal errors caused 
    by `DelaunayTriangulation.triangulate`, file an issue at 
    [DelaunayTriangulation.jl](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/issues/new).

## References

* Cheng et al. 2012. [Delaunay Mesh Generation]
  (https://people.eecs.berkeley.edu/~jrs/meshbook.html)
"""
struct DelaunayTriangulation{RNG<:AbstractRNG} <: BoundaryDiscretizationMethod
  rng::RNG
end

DelaunayTriangulation(rng=Random.default_rng()) = DelaunayTriangulation(rng)

function discretizewithin(ring::Ring{2}, method::DelaunayTriangulation)
  points = vertices(ring)
  coords = map(p -> ustrip.(to(p)), points)
  bnodes = [1:nvertices(ring); 1]
  triang = triangulate(coords, boundary_nodes=bnodes, rng=method.rng)
  connec = connect.(each_solid_triangle(triang))
  SimpleMesh(points, connec)
end
