using BenchmarkTools
using Meshes
using Unitful
using LinearAlgebra

# auxiliary variables
point = Point(0, 0)
points = rand(Point, 100)
sphere = Sphere((0, 0), 1)
poly1 = PolyArea([sphere(t) for t in 0.1:0.1:1.0])
poly2 = PolyArea([sphere(t) + Vec(1, 0) for t in 0.1:0.1:1.0])
ring1 = Ring([sphere(t) for t in range(0.1, 1.0, length=5)])
ring2 = Ring([sphere(t) for t in range(0.1, 1.0, length=1500)])
grid = CartesianGrid(100, 100)
mesh = discretize(Sphere((0, 0, 0), 1))
ray = Ray((-1, -1, -1), (0, 0, 1))
tri = Triangle((0, 0, 0), (1, 0, 0), (0, 1, 0))
tet = Tetrahedron((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
poly = PolyArea((0, 0), (1, 0), (1, 1), (0.5, 2), (0, 1))
search1 = KNearestSearch(grid, 20)
search2 = BallSearch(grid, MetricBall(30))
search3 = KBallSearch(grid, 20, MetricBall(30))
search4 = BoundingBoxSearch(grid)

# initialize benchmark suite
const SUITE = BenchmarkGroup()

# ---------
# CLIPPING
# ---------

SUITE["clipping"] = BenchmarkGroup()

SUITE["clipping"]["SutherlandHodgman"] = @benchmarkable clip($poly1, $poly2, SutherlandHodgmanClipping())

# ---------------
# DISCRETIZATION
# ---------------

SUITE["discretization"] = BenchmarkGroup()

SUITE["discretization"]["simplexify"] = @benchmarkable simplexify($mesh)

# ---------
# TOPOLOGY
# ---------

SUITE["topology"] = BenchmarkGroup()

SUITE["topology"]["half-edge"] = @benchmarkable convert(HalfEdgeTopology, topology($mesh))

# ----------------
# NEIGHBOR SEARCH
# ----------------

SUITE["neighborsearch"] = BenchmarkGroup()

SUITE["neighborsearch"]["knn"] = @benchmarkable search($point, $search1)
SUITE["neighborsearch"]["ball"] = @benchmarkable search($point, $search2)
SUITE["neighborsearch"]["knnball"] = @benchmarkable search($point, $search3)
SUITE["neighborsearch"]["bbox-point"] = @benchmarkable search($point, $search4)
SUITE["neighborsearch"]["bbox-poly"] = @benchmarkable search($poly, $search4)

# --------
# WINDING
# --------

SUITE["winding"] = BenchmarkGroup()

SUITE["winding"]["mesh"] = @benchmarkable winding($points, $mesh)

# -------
# SIDEOF
# -------

SUITE["sideof"] = BenchmarkGroup()

SUITE["sideof"]["ring"]["small"] = @benchmarkable sideof($point, $ring1)
SUITE["sideof"]["ring"]["large"] = @benchmarkable sideof($point, $ring2)

# -----------
# INTERSECTS
# -----------

SUITE["intersects"] = BenchmarkGroup()

SUITE["intersects"]["triangle-tetrahedron"] = @benchmarkable intersects($tri, $tet)

# -------------
# INTERSECTION
# -------------

SUITE["intersection"] = BenchmarkGroup()

SUITE["intersection"]["ray-triangle"] = @benchmarkable intersection($ray, $tri)

# -------------
# DIFFERENTIATION
# -------------

SUITE["differentiation"] = BenchmarkGroup()

SUITE["differentiation"]["ray"] = @benchmarkable jacobian($ray, (0.2,))
SUITE["differentiation"]["triangle"] = @benchmarkable jacobian($tri, (0.2, 0.3))
SUITE["differentiation"]["tetrahedron"] = @benchmarkable jacobian($tet, (0.2, 0.3, 0.4))

# -------------
# INTEGRATION
# -------------

SUITE["integration"] = BenchmarkGroup()

SUITE["integration"]["ray"] = @benchmarkable integral($ray) do p
  r = ustrip(norm(to(p)))
  exp(-r^2) * u"A"
end
SUITE["integration"]["triangle"] = @benchmarkable integral(p -> 1, $tri)
SUITE["integration"]["tetrahedron"] = @benchmarkable integral(p -> 1, $tet)

# ---------
# CENTROID
# ---------

SUITE["centroid"] = BenchmarkGroup()

SUITE["centroid"]["polyarea"] = @benchmarkable centroid($poly)
