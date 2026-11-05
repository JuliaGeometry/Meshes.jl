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
  DelaunayRefinement(rng, asdeg(minangle), isnothing(maxarea) ? nothing : _asarea(maxarea))

function refine(mesh, method::DelaunayRefinement)
  assertion(embeddim(mesh) == 2 && manifold(mesh) <: 𝔼, "DelaunayRefinement only defined for 2D meshes")

  # triangles without repeated vertices, as merged meshes have them
  tmesh = simplexify(mesh |> Repair(0))
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

_asarea(x::Quantity{<:Any,u"𝐋^2"}) = x
_asarea(x::Number) = x * u"m^2"
_asarea(::Quantity) = throw(ArgumentError("invalid area unit"))

# closed curves along the boundary of a triangle mesh, counter-clockwise
# around the domain and clockwise around holes, the outer curve first
function _boundarycurves(mesh)
  points = vertices(mesh)

  # edges used by one triangle only, oriented with the triangle
  nuses = Dict{NTuple{2,Int},Int}()
  edges = NTuple{2,Int}[]
  for elem in elements(topology(mesh))
    i, j, k = indices(elem)
    isnegative(signarea(points[i], points[j], points[k])) && ((j, k) = (k, j))
    for e in ((i, j), (j, k), (k, i))
      push!(edges, e)
      nuses[minmax(e...)] = get(nuses, minmax(e...), 0) + 1
    end
  end
  succ = Dict{Int,Int}()
  for (i, j) in edges
    nuses[minmax(i, j)] == 1 || continue
    assertion(!haskey(succ, i), "DelaunayRefinement needs a boundary that does not touch itself")
    succ[i] = j
  end

  # chain the boundary edges into closed curves
  curves = Vector{Int}[]
  seen = Set{Int}()
  for s in sort!(collect(keys(succ)))
    s ∈ seen && continue
    curve = [s]
    v = succ[s]
    while v ≠ s
      push!(curve, v)
      push!(seen, v)
      v = succ[v]
    end
    push!(seen, s)
    push!(curves, [curve; s])
  end

  # the outer curve has the largest area
  sarea(c) = sum(i -> signarea(points[c[1]], points[c[i]], points[c[i + 1]]), 2:(length(c) - 2))
  sort!(curves, by=sarea, rev=true)
end
