# ------------------------------------------------------------------
# Licensed under the MIT License. See LICENSE in the project root.
# ------------------------------------------------------------------

"""
    DelaunayRefinement([rng]; minangle=30, maxarea=nothing)

Delaunay refinement of 2D meshes into triangles.
Points are added until the triangles have angles of at least `minangle` (in degrees by default).
If `maxarea` is given, the triangles also have an area of at most `maxarea` (in square meters by default).
Corners of the domain with a smaller angle than `minangle` keep their angle.
Optionally, specify the random number generator `rng`.

## References

* Cheng et al. 2012. [Delaunay Mesh Generation](https://people.eecs.berkeley.edu/~jrs/meshbook.html)

### Notes

Wraps the `refine!` function of DelaunayTriangulation.jl.
With `minangle` above about 33.9 degrees, the refinement may not converge and can create millions of triangles.
For any internal errors, file an issue at [DelaunayTriangulation.jl](https://github.com/JuliaGeometry/DelaunayTriangulation.jl/issues/new)
"""
struct DelaunayRefinement{RNG<:AbstractRNG,ℳ<:Deg,𝒜} <: RefinementMethod
  rng::RNG
  minangle::ℳ
  maxarea::𝒜
end

DelaunayRefinement(rng=Random.default_rng(); minangle=30, maxarea=nothing) =
  DelaunayRefinement(rng, asdeg(minangle), isnothing(maxarea) ? nothing : asarea(maxarea))

function refine(mesh, method::DelaunayRefinement)
  assertion(embeddim(mesh) == 2 && manifold(mesh) <: 𝔼, "DelaunayRefinement only defined for 2D meshes")

  tmesh = simplexify(mesh)
  points = vertices(tmesh)
  ℒ = lentype(tmesh)
  u = unit(ℒ)
  coords = map(p -> Tuple(ustrip.(u, to(p))), points)

  # Delaunay triangulation of the domain, refined with the constraints
  curves = _boundarycurves(tmesh)
  triang = triangulate(coords, boundary_nodes=[[curve] for curve in curves], rng=method.rng)
  min_angle = ustrip(u"°", method.minangle)
  max_area = isnothing(method.maxarea) ? typemax(numtype(ℒ)) : ustrip(u^2, method.maxarea)
  refine!(triang; min_angle, max_area, rng=method.rng)

  # original points are kept as they are, new points are added
  inds = sort!(collect(each_solid_vertex(triang)))
  newind = Dict(v => k for (k, v) in enumerate(inds))
  newpoints = map(inds) do v
    v ≤ length(points) ? points[v] : withcrs(tmesh, getxy(get_point(triang, v)) .* u)
  end
  connec = map(t -> connect(map(v -> newind[v], t)), each_solid_triangle(triang))

  SimpleMesh(newpoints, connec)
end

# closed curves of vertex indices along the boundary of a mesh, counter-clockwise
# around the domain and clockwise around holes, the outer curve first
function _boundarycurves(mesh)
  bound = boundary(mesh)
  brings = bound isa Multi ? parent(bound) : [bound]

  # the outer curve encloses the holes, so it has the largest area
  curves = sort(brings, by=r -> measure(Pgon(r)), rev=true)
  outer = first(curves)

  # a curve outside of the outer one bounds another piece of the mesh
  insideouter(r) = all(v -> sideof(v, outer) ≠ OUT, eachvertex(r))
  assertion(all(insideouter, curves[2:end]), "DelaunayRefinement needs a mesh with a single outer boundary")

  index = Dict(v => i for (i, v) in enumerate(vertices(mesh)))
  map(enumerate(curves)) do (i, r)
    ring = orientation(r) == (isone(i) ? CCW : CW) ? r : reverse(r)
    inds = [index[v] for v in eachvertex(ring)]
    [inds; first(inds)]
  end
end
