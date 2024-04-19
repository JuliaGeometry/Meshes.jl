@testset "Sampling" begin
  @testset "UniformSampling" begin
    rng = StableRNG(123)
    d = CartesianGrid{T}(100, 100)
    s = sample(rng, d, UniformSampling(100))
    μ = mean(coordinates.([centroid(s, i) for i in 1:nelements(s)]))
    @test nelements(s) == 100
    @test isapprox(μ, T[50.0, 50.0], atol=T(10))

    # availability of option ordered
    s = sample(rng, d, UniformSampling(100, ordered=true))
    μ = mean(coordinates.([centroid(s, i) for i in 1:nelements(s)]))
    @test nelements(s) == 100
    @test isapprox(μ, T[50.0, 50.0], atol=T(10))
  end

  @testset "WeightedSampling" begin
    # uniform weights => uniform sampler
    rng = StableRNG(123)
    d = CartesianGrid{T}(100, 100)
    s = sample(rng, d, WeightedSampling(100))
    μ = mean(coordinates.([centroid(s, i) for i in 1:nelements(s)]))
    @test nelements(s) == 100
    @test isapprox(μ, T[50.0, 50.0], atol=T(10))

    # availability of option ordered
    s = sample(rng, d, WeightedSampling(100, ordered=true))
    μ = mean(coordinates.([centroid(s, i) for i in 1:nelements(s)]))
    @test nelements(s) == 100
    @test isapprox(μ, T[50.0, 50.0], atol=T(10))

    # utility method
    s = sample(rng, d, 100, ordered=true)
    μ = mean(coordinates.([centroid(s, i) for i in 1:nelements(s)]))
    @test nelements(s) == 100
    @test isapprox(μ, T[50.0, 50.0], atol=T(10))
    s = sample(rng, d, 100, fill(1, 10000), ordered=true)
    μ = mean(coordinates.([centroid(s, i) for i in 1:nelements(s)]))
    @test nelements(s) == 100
    @test isapprox(μ, T[50.0, 50.0], atol=T(10))
  end

  @testset "BallSampling" begin
    d = CartesianGrid{T}(100, 100)
    s = sample(d, BallSampling(T(10)))
    n = nelements(s)
    x = coordinates(centroid(s, 1))
    y = coordinates(centroid(s, 17))
    @test n < 100
    @test sqrt(sum((x - y) .^ 2)) ≥ T(10)

    d = CartesianGrid{T}(100, 100)
    s = sample(d, BallSampling(T(20)))
    n = nelements(s)
    x = coordinates(centroid(s, 1))
    y = coordinates(centroid(s, 17))
    @test n < 50
    @test sqrt(sum((x - y) .^ 2)) ≥ T(20)
  end

  @testset "BlockSampling" begin
    g = CartesianGrid{T}(100, 100)
    s = sample(g, BlockSampling(T(10)))
    @test nelements(s) == 100
    x = coordinates.(centroid.(s))
    D = pairwise(Euclidean(), x)
    d = [D[i, j] for i in 1:length(x) for j in 1:(i - 1)]
    @test all(≥(T(10)), d)
  end

  @testset "RegularSampling" begin
    b = Box(P2(0, 0), P2(2, 2))
    ps = sample(b, RegularSampling(3))
    @test collect(ps) == P2[(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (0, 2), (1, 2), (2, 2)]
    ps = sample(b, RegularSampling(2, 3))
    @test collect(ps) == P2[(0, 0), (2, 0), (0, 1), (2, 1), (0, 2), (2, 2)]

    b = BezierCurve([P2(0, 0), P2(1, 0), P2(1, 1)])
    ps = sample(b, RegularSampling(4))
    ts = P2[(0.0, 0.0), (0.5555555555555556, 0.1111111111111111), (0.8888888888888888, 0.4444444444444444), (1.0, 1.0)]
    for (p, t) in zip(ps, ts)
      @test p ≈ t
    end

    s = Sphere(P2(0, 0), T(2))
    ps = sample(s, RegularSampling(4))
    ts = P2[(2, 0), (0, 2), (-2, 0), (0, -2)]
    for (p, t) in zip(ps, ts)
      @test p ≈ t
    end

    s = Sphere(P3(0, 0, 0), T(2))
    ps = sample(s, RegularSampling(2, 2))
    ts = P3[
      (1.7320508075688772, 0.0, 1.0),
      (1.7320508075688772, 0.0, -1.0),
      (-1.7320508075688772, 0.0, 1.0),
      (-1.7320508075688772, 0.0, -1.0),
      (0.0, 0.0, 2.0),
      (0.0, 0.0, -2.0)
    ]
    for (p, t) in zip(ps, ts)
      @test p ≈ t
    end

    e = Ellipsoid((T(3), T(2), T(1)), P3(1, 1, 1), RotZYX(T(π / 4), T(π / 4), T(π / 4)))
    ps = sample(e, RegularSampling(2, 2))
    ts = P3[
      (2.725814800973295, 2.225814800973295, -0.5871173070873834),
      (1.872261410380021, 2.372261410380021, -1.0871173070873832),
      (0.12773858961997864, -0.37226141038002103, 3.0871173070873836),
      (-0.725814800973295, -0.22581480097329454, 2.587117307087383),
      (1.8535533905932737, 0.8535533905932737, 1.5),
      (0.14644660940672627, 1.1464466094067263, 0.4999999999999999)
    ]
    for (p, t) in zip(ps, ts)
      @test p ≈ t
    end

    b = Ball(P2(0, 0), T(2))
    ps = sample(b, RegularSampling(3, 4))
    @test all(∈(b), ps)
    ts = P2[
      (0.6666666666666666, 0.0),
      (1.3333333333333333, 0.0),
      (2.0, 0.0),
      (0.0, 0.6666666666666666),
      (0.0, 1.3333333333333333),
      (0.0, 2.0),
      (-0.6666666666666666, 0.0),
      (-1.3333333333333333, 0.0),
      (-2.0, 0.0),
      (0.0, -0.6666666666666666),
      (0.0, -1.3333333333333333),
      (0.0, -2.0),
      (0.0, 0.0)
    ]
    for (p, t) in zip(ps, ts)
      @test p ≈ t
    end

    b = Ball(P2(10, 10), T(2))
    ps = sample(b, RegularSampling(4, 3))
    @test all(∈(b), ps)
    ts = P2[
      (10.5, 10.0),
      (11.0, 10.0),
      (11.5, 10.0),
      (12.0, 10.0),
      (9.75, 10.433012701892219),
      (9.5, 10.86602540378444),
      (9.25, 11.299038105676658),
      (9.0, 11.732050807568877),
      (9.75, 9.566987298107781),
      (9.5, 9.13397459621556),
      (9.25, 8.700961894323342),
      (9.0, 8.267949192431121),
      (10.0, 10.0)
    ]
    for (p, t) in zip(ps, ts)
      @test p ≈ t
    end

    b = Ball(P3(0, 0, 0), T(2))
    ps = sample(b, RegularSampling(3, 2, 3))
    @test all(∈(b), ps)
    ts = P3[
      (0.5773502691896257, 0.0, 0.3333333333333333),
      (1.1547005383792515, 0.0, 0.6666666666666666),
      (1.7320508075688772, 0.0, 1.0),
      (0.5773502691896256, 0.0, -0.3333333333333335),
      (1.1547005383792512, 0.0, -0.666666666666667),
      (1.732050807568877, 0.0, -1.0000000000000004),
      (-0.288675134594813, 0.4999999999999999, 0.3333333333333333),
      (-0.577350269189626, 0.9999999999999998, 0.6666666666666666),
      (-0.8660254037844389, 1.4999999999999996, 1.0),
      (-0.2886751345948129, 0.4999999999999998, -0.3333333333333335),
      (-0.5773502691896258, 0.9999999999999996, -0.666666666666667),
      (-0.8660254037844388, 1.4999999999999993, -1.0000000000000004),
      (-0.28867513459481264, -0.5000000000000001, 0.3333333333333333),
      (-0.5773502691896253, -1.0000000000000002, 0.6666666666666666),
      (-0.8660254037844379, -1.5000000000000004, 1.0),
      (-0.2886751345948126, -0.5, -0.3333333333333335),
      (-0.5773502691896252, -1.0, -0.666666666666667),
      (-0.8660254037844378, -1.5000000000000002, -1.0000000000000004),
      (0.0, 0.0, 0.0)
    ]
    for (p, t) in zip(ps, ts)
      @test p ≈ t
    end

    b = Ball(P3(10, 10, 10), T(2))
    ps = sample(b, RegularSampling(3, 2, 3))
    @test all(∈(b), ps)

    # cylinder surface with parallel planes
    c = CylinderSurface(Plane(P3(0, 0, 0), V3(0, 0, 1)), Plane(P3(0, 0, 1), V3(0, 0, 1)), T(1))
    ps = sample(c, RegularSampling(20, 10))
    cs = coordinates.(ps)
    xs = getindex.(cs, 1)
    ys = getindex.(cs, 2)
    zs = getindex.(cs, 3)
    @test length(cs) == 200 + 2
    @test all(T(-1) ≤ x ≤ T(1) for x in xs)
    @test all(T(-1) ≤ y ≤ T(1) for y in ys)
    @test all(T(0) ≤ z ≤ T(1) for z in zs)

    # cylinder surface with parallel shifted planes
    c = CylinderSurface(Plane(P3(0, 0, 0), V3(0, 0, 1)), Plane(P3(1, 1, 1), V3(0, 0, 1)), T(1))
    ps = sample(c, RegularSampling(20, 10))
    cs = coordinates.(ps)
    xs = getindex.(cs, 1)
    ys = getindex.(cs, 2)
    zs = getindex.(cs, 3)
    @test length(cs) == 200 + 2

    # cylinder surface with non-parallel planes
    c = CylinderSurface(Plane(P3(0, 0, 0), V3(1, 0, 1)), Plane(P3(1, 1, 1), V3(0, 1, 1)), T(1))
    ps = sample(c, RegularSampling(20, 10))
    cs = coordinates.(ps)
    @test length(cs) == 200 + 2

    s = Segment(P2(0, 0), P2(1, 1))
    ps = sample(s, RegularSampling(2))
    @test collect(ps) == P2[(0, 0), (1, 1)]
    ps = sample(s, RegularSampling(3))
    @test collect(ps) == P2[(0, 0), (0.5, 0.5), (1, 1)]

    q = Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
    ps = sample(q, RegularSampling(2, 2))
    @test collect(ps) == P2[(0, 0), (1, 0), (0, 1), (1, 1)]
    ps = sample(q, RegularSampling(3, 3))
    @test collect(ps) == P2[(0, 0), (0.5, 0), (1, 0), (0, 0.5), (0.5, 0.5), (1, 0.5), (0, 1), (0.5, 1), (1, 1)]

    h =
      Hexahedron(P3(0, 0, 0), P3(1, 0, 0), P3(1, 1, 0), P3(0, 1, 0), P3(0, 0, 1), P3(1, 0, 1), P3(1, 1, 1), P3(0, 1, 1))
    ps = sample(h, RegularSampling(2, 2, 2))
    @test collect(ps) == P3[(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 1)]
    ps = sample(h, RegularSampling(3, 2, 2))
    @test collect(ps) == P3[
      (0, 0, 0),
      (0.5, 0, 0),
      (1, 0, 0),
      (0, 1, 0),
      (0.5, 1, 0),
      (1, 1, 0),
      (0, 0, 1),
      (0.5, 0, 1),
      (1, 0, 1),
      (0, 1, 1),
      (0.5, 1, 1),
      (1, 1, 1)
    ]

    torus = Torus(P3(0, 0, 0), V3(1, 0, 0), T(2), T(1))
    ps = sample(torus, RegularSampling(3, 3))
    ts = P3[
      (0, 0, -3),
      (-sqrt(3) / 2, 0, -1.5),
      (sqrt(3) / 2, 0, -1.5),
      (0, 3sqrt(3) / 2, 1.5),
      (-sqrt(3) / 2, 3sqrt(3) / 4, 0.75),
      (sqrt(3) / 2, 3sqrt(3) / 4, 0.75),
      (0, -3sqrt(3) / 2, 1.5),
      (-sqrt(3) / 2, -3sqrt(3) / 4, 0.75),
      (sqrt(3) / 2, -3sqrt(3) / 4, 0.75)
    ]
    for (p, t) in zip(ps, ts)
      @test p ≈ t
    end

    grid = CartesianGrid{T}(10, 10)
    points = sample(grid, RegularSampling(100, 200))
    @test length(collect(points)) == 20000
  end

  @testset "HomogeneousSampling" begin
    s = Segment(P2(0, 0), P2(1, 0))
    ps = sample(s, HomogeneousSampling(100))
    @test first(ps) isa P2
    @test all(0 ≤ coords[1] ≤ 1 for coords in coordinates.(ps))
    @test all(coords[2] == 0 for coords in coordinates.(ps))

    s = Segment(P2(0, 0), P2(0, 1))
    ps = sample(s, HomogeneousSampling(100))
    @test first(ps) isa P2
    @test all(coords[1] == 0 for coords in coordinates.(ps))
    @test all(0 ≤ coords[2] ≤ 1 for coords in coordinates.(ps))

    s = Segment(P2(0, 0), P2(1, 1))
    ps = sample(s, HomogeneousSampling(100))
    @test first(ps) isa P2
    @test all(0 ≤ coords[1] == coords[2] ≤ 1 for coords in coordinates.(ps))

    c = Rope(P2(0, 0), P2(1, 0), P2(0, 1), P2(1, 1))
    ps = sample(c, HomogeneousSampling(100))
    @test first(ps) isa P2
    @test all(coords[1] + coords[2] == 1 || (0 ≤ coords[1] ≤ 1 && coords[2] ∈ [0, 1]) for coords in coordinates.(ps))

    t = Triangle(P2(0, 0), P2(1, 0), P2(0, 1))
    ps = sample(t, HomogeneousSampling(100))
    @test first(ps) isa P2
    @test all(∈(t), ps)

    q = Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
    ps = sample(q, HomogeneousSampling(100))
    @test first(ps) isa P2
    @test all(∈(q), ps)

    b = Ball(P2(10, 10), T(3))
    ps = sample(b, HomogeneousSampling(100))
    @test first(ps) isa P2
    @test all(∈(b), ps)

    b = Ball(P3(10, 10, 10), T(10))
    ps = sample(b, HomogeneousSampling(100))
    @test first(ps) isa P3
    @test all(∈(b), ps)

    poly1 = PolyArea(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
    poly2 = PolyArea(P2[(1, 1), (2, 1), (2, 2), (1, 2)])
    multi = Multi([poly1, poly2])
    ps = sample(multi, HomogeneousSampling(100))
    @test all(p -> (P2(0, 0) ⪯ p ⪯ P2(1, 1)) || (P2(1, 1) ⪯ p ⪯ P2(2, 2)), ps)

    points = P2[(0, 0), (1, 0), (0, 1), (1, 1), (0.25, 0.5), (0.75, 0.5)]
    connec = connect.([(3, 1, 5), (4, 6, 2), (1, 2, 6, 5), (5, 6, 4, 3)])
    mesh = SimpleMesh(points, connec)
    ps = sample(mesh, HomogeneousSampling(400))
    @test first(ps) isa P2
    @test all(∈(mesh), ps)
    ps = sample(mesh, HomogeneousSampling(400, 1:nelements(mesh)))
    @test first(ps) isa P2
    @test all(∈(mesh), ps)
  end

  @testset "MinDistanceSampling" begin
    poly1 = PolyArea(P2[(0, 0), (1, 0), (1, 1), (0, 1)])
    poly2 = PolyArea(P2[(1, 1), (2, 1), (2, 2), (1, 2)])
    multi = Multi([poly1, poly2])
    ps = sample(multi, MinDistanceSampling(0.1))
    @test all(p -> (P2(0, 0) ⪯ p ⪯ P2(1, 1)) || (P2(1, 1) ⪯ p ⪯ P2(2, 2)), ps)

    points = P2[(0, 0), (1, 0), (0, 1), (1, 1), (0.25, 0.5), (0.75, 0.5)]
    connec = connect.([(3, 1, 5), (4, 6, 2), (1, 2, 6, 5), (5, 6, 4, 3)])
    mesh = SimpleMesh(points, connec)
    ps = sample(mesh, MinDistanceSampling(0.2))
    n = length(ps)
    @test first(ps) isa P2
    @test all(∈(mesh), ps)
    @test all(norm(ps[i] - ps[j]) ≥ 0.2 for i in 1:n for j in (i + 1):n)

    # geometries with almost zero measure
    # can still be sampled (at least one point)
    poly = PolyArea(P2[(-44.20065308, -21.12284851), (-44.20324135, -21.122799875), (-44.20582962, -21.12275124)])
    ps = sample(poly, MinDistanceSampling(3.2423333333753135e-5))
    @test length(ps) > 0
  end

  @testset "RNGs" begin
    dom = CartesianGrid{T}(100, 100)
    for method in [UniformSampling(100), WeightedSampling(100), BallSampling(T(10))]
      rng = StableRNG(2021)
      s1 = sample(rng, dom, method)
      rng = StableRNG(2021)
      s2 = sample(rng, dom, method)
      @test collect(s1) == collect(s2)
    end

    # cannot test some sampling methods with T = Float32
    # because of https://github.com/JuliaStats/StatsBase.jl/issues/695
    if T == Float64
      for method in [HomogeneousSampling(100), MinDistanceSampling(T(5))]
        rng = StableRNG(2021)
        s1 = sample(rng, dom, method)
        rng = StableRNG(2021)
        s2 = sample(rng, dom, method)
        @test collect(s1) == collect(s2)
      end
    end

    method = RegularSampling(10)
    for geom in [
      Box(P2(0, 0), P2(2, 2))
      Sphere(P2(0, 0), T(2))
      Ball(P2(0, 0), T(2))
      Segment(P2(0, 0), P2(1, 1))
      Quadrangle(P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1))
      Hexahedron(
        P3(0, 0, 0),
        P3(1, 0, 0),
        P3(1, 1, 0),
        P3(0, 1, 0),
        P3(0, 0, 1),
        P3(1, 0, 1),
        P3(1, 1, 1),
        P3(0, 1, 1)
      )
    ]
      rng = StableRNG(2021)
      s1 = sample(rng, geom, method)
      rng = StableRNG(2021)
      s2 = sample(rng, geom, method)
      @test collect(s1) == collect(s2)
    end
  end
end
