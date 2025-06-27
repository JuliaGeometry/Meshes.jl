using BenchmarkTools
using Meshes

# auxiliary variables for tests
pₒ = Point(0, 0)
ps = rand(Point, 100)
s = Sphere((0, 0), 1)
r₁ = Ring([s(t) for t in 0.1:0.1:1.0])
r₂ = Ring([s(t) + Vec(1, 0) for t in 0.1:0.1:1.0])
rₛ = Ring([s(t) for t in range(0.1, 1.0, length=5)])
rₗ = Ring([s(t) for t in range(0.1, 1.0, length=1500)])
m = discretize(Sphere((0, 0, 0), 1))

# initialize benchmark suite
const SUITE = BenchmarkGroup()

# ---------
# CLIPPING
# ---------

SUITE["clipping"] = BenchmarkGroup()

SUITE["clipping"]["SutherlandHodgman"] = @benchmarkable clip($r₁, $r₂, SutherlandHodgmanClipping())

# ---------------
# DISCRETIZATION
# ---------------

SUITE["discretization"] = BenchmarkGroup()

SUITE["discretization"]["simplexify"] = @benchmarkable simplexify($m)

# --------
# WINDING
# --------

SUITE["winding"] = BenchmarkGroup()

SUITE["winding"]["mesh"] = @benchmarkable winding($ps, $m)

# -------
# SIDEOF
# -------

SUITE["sideof"] = BenchmarkGroup()

SUITE["sideof"]["ring"]["small"] = @benchmarkable sideof($pₒ, $rₛ)
SUITE["sideof"]["ring"]["large"] = @benchmarkable sideof($pₒ, $rₗ)
