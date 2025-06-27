using BenchmarkTools
using Meshes

const SUITE = BenchmarkGroup()

# ---------
# CLIPPING
# ---------

SUITE["clipping"] = BenchmarkGroup()

s1 = Sphere((0, 0), 1)
s2 = Sphere((1, 0), 1)
r1 = Ring([s1(t) for t in 0.1:0.1:1.0])
r2 = Ring([s2(t) for t in 0.1:0.1:1.0])

SUITE["clipping"]["SutherlandHodgman"] = @benchmarkable clip($r1, $r2, SutherlandHodgmanClipping())

# ---------------
# DISCRETIZATION
# ---------------

SUITE["discretization"] = BenchmarkGroup()

m = discretize(Sphere(0, 0, 0), 1)

SUITE["discretization"]["simplexify"] = @benchmarkable simplexify($m)
