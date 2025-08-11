using BenchmarkTools
using Meshes
using Random

# auxiliary variables
point1 = Point(0, 0)
points = rand(Point, 100)
sphere = Sphere((0, 0), 1)
ring1 = Ring([sphere(t) for t in 0.1:0.1:1.0])
ring2 = Ring([sphere(t) + Vec(1, 0) for t in 0.1:0.1:1.0])
ring3 = Ring([sphere(t) for t in range(0.1, 1.0, length=5)])
ring4 = Ring([sphere(t) for t in range(0.1, 1.0, length=1500)])
mesh = discretize(Sphere((0, 0, 0), 1))
bigmesh = refine(refine(mesh, TriRefinement()), TriSubdivision())
bigmesh_connec = collect(elements(topology(bigmesh)))
ray = Ray((-1, -1, -1), (0, 0, 1))
triangle = Triangle((0, 0, 0), (1, 0, 0), (0, 1, 0))

# initialize benchmark suite
const SUITE = BenchmarkGroup()

# ---------
# CLIPPING
# ---------

SUITE["clipping"] = BenchmarkGroup()

SUITE["clipping"]["SutherlandHodgman"] = @benchmarkable clip($ring1, $ring2, SutherlandHodgmanClipping())

# ---------------
# DISCRETIZATION
# ---------------

SUITE["discretization"] = BenchmarkGroup()

SUITE["discretization"]["simplexify"] = @benchmarkable simplexify($mesh)

# ---------
# TOPOLOGY
# ---------

SUITE["topology"] = BenchmarkGroup()

SUITE["topology"]["half-edge"] = BenchmarkGroup()

SUITE["topology"]["half-edge"]["adjsortperm-defaultorder"] =
  @benchmarkable Meshes.adjsortperm(connec) setup = (connec = bigmesh_connec)
SUITE["topology"]["half-edge"]["adjsortperm-shuffled"] =
  @benchmarkable Meshes.adjsortperm(connec) setup = (connec = shuffle(bigmesh_connec))
SUITE["topology"]["half-edge"]["adjsortperm-presorted"] =
  @benchmarkable Meshes.adjsortperm(connec) setup = (connec = bigmesh_connec[Meshes.adjsortperm(bigmesh_connec)])
SUITE["topology"]["half-edge"]["adjsortperm-reverse-sorted"] = @benchmarkable Meshes.adjsortperm(connec) setup = begin
  sp = Meshes.adjsortperm(bigmesh_connec)
  connec = bigmesh_connec[sp[[1; end:-1:2]]]
end

SUITE["topology"]["half-edge"]["constructor-defaultorder", "sort=false"] =
  @benchmarkable HalfEdgeTopology(connec; sort=false) setup = (connec = bigmesh_connec)
SUITE["topology"]["half-edge"]["constructor-shuffled", "sort=false"] =
  @benchmarkable HalfEdgeTopology(connec; sort=false) setup = (connec = shuffle(bigmesh_connec))
SUITE["topology"]["half-edge"]["constructor-presorted", "sort=false"] =
  @benchmarkable HalfEdgeTopology(connec; sort=false) setup =
    (connec = bigmesh_connec[Meshes.adjsortperm(bigmesh_connec)])
SUITE["topology"]["half-edge"]["constructor-reverse-sorted", "sort=false"] =
  @benchmarkable HalfEdgeTopology(connec; sort=false) setup = begin
    sp = Meshes.adjsortperm(bigmesh_connec)
    connec = bigmesh_connec[sp[[1; end:-1:2]]]
  end

SUITE["topology"]["half-edge"]["constructor-defaultorder", "sort=true"] =
  @benchmarkable HalfEdgeTopology(connec; sort=true) setup = (connec = bigmesh_connec)
SUITE["topology"]["half-edge"]["constructor-shuffled", "sort=true"] =
  @benchmarkable HalfEdgeTopology(connec; sort=true) setup = (connec = shuffle(bigmesh_connec))
SUITE["topology"]["half-edge"]["constructor-presorted", "sort=true"] =
  @benchmarkable HalfEdgeTopology(connec; sort=true) setup =
    (connec = bigmesh_connec[Meshes.adjsortperm(bigmesh_connec)])
SUITE["topology"]["half-edge"]["constructor-reverse-sorted", "sort=true"] =
  @benchmarkable HalfEdgeTopology(reverse(connec); sort=true) setup = begin
    sp = Meshes.adjsortperm(bigmesh_connec)
    connec = bigmesh_connec[sp[[1; end:-1:2]]]
  end

# --------
# WINDING
# --------

SUITE["winding"] = BenchmarkGroup()

SUITE["winding"]["mesh"] = @benchmarkable winding($points, $mesh)

# -------
# SIDEOF
# -------

SUITE["sideof"] = BenchmarkGroup()

SUITE["sideof"]["ring"]["small"] = @benchmarkable sideof($point1, $ring3)
SUITE["sideof"]["ring"]["large"] = @benchmarkable sideof($point1, $ring4)

# -------------
# INTERSECTION
# -------------

SUITE["intersection"] = BenchmarkGroup()

SUITE["intersection"]["ray-triangle"] = @benchmarkable intersection($ray, $triangle)
