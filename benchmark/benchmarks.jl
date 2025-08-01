using BenchmarkTools
using Meshes

# auxiliary variables
point1 = Point(0, 0)
points = rand(Point, 100)
sphere = Sphere((0, 0), 1)
ring1 = Ring([sphere(t) for t in 0.1:0.1:1.0])
ring2 = Ring([sphere(t) + Vec(1, 0) for t in 0.1:0.1:1.0])
ring3 = Ring([sphere(t) for t in range(0.1, 1.0, length=5)])
ring4 = Ring([sphere(t) for t in range(0.1, 1.0, length=1500)])
mesh = discretize(Sphere((0, 0, 0), 1))
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

# --------
# TOPOLOGIES
# --------

SUITE["topologies"] = BenchmarkGroup()

SUITE["topologies"]["half-edge"] = @benchmarkable convert(HalfEdgeTopology, topology($mesh))

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
