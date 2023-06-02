@testset "Sampling" begin
  @testset "UniformSampling" begin
    Random.seed!(2021)
    d = CartesianGrid{T}(100, 100)
    s = sample(d, UniformSampling(100))
    μ = mean(coordinates.([centroid(s, i) for i in 1:nelements(s)]))
    @test nelements(s) == 100
    @test isapprox(μ, T[50.0, 50.0], atol=T(10))

    # availability of option ordered
    s = sample(d, UniformSampling(100, ordered=true))
    μ = mean(coordinates.([centroid(s, i) for i in 1:nelements(s)]))
    @test nelements(s) == 100
    @test isapprox(μ, T[50.0, 50.0], atol=T(10))
  end

  @testset "WeightedSampling" begin
    # uniform weights => uniform sampler
    Random.seed!(2020)
    d = CartesianGrid{T}(100, 100)
    s = sample(d, WeightedSampling(100))
    μ = mean(coordinates.([centroid(s, i) for i in 1:nelements(s)]))
    @test nelements(s) == 100
    @test isapprox(μ, T[50.0, 50.0], atol=T(10))

    # availability of option ordered
    s = sample(d, WeightedSampling(100, ordered=true))
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
    # fix import conflict with Plots
    BezierCurve = Meshes.BezierCurve

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
      (-1.7320508075688772, 0.0, -1.0)
    ]
    for (p, t) in zip(ps, ts)
      @test p ≈ t
    end

    b = Ball(P2(0, 0), T(2))
    ps = sample(b, RegularSampling(4, 3))
    ts = P2[
      (1.0, 0.0),
      (0.0, 1.0),
      (-1.0, 0.0),
      (0.0, -1.0),
      (1.5, 0.0),
      (0.0, 1.5),
      (-1.5, 0.0),
      (0.0, -1.5),
      (2.0, 0.0),
      (0.0, 2.0),
      (-2.0, 0.0),
      (0.0, -2.0)
    ]
    for (p, t) in zip(ps, ts)
      @test p ≈ t
    end

    b = Ball(P2(10, 10), T(2))
    ps = sample(b, RegularSampling(4, 3))
    @test all(∈(b), ps)
    ts = P2[
      (11.0, 10.0),
      (10.0, 11.0),
      (9.0, 10.0),
      (10.0, 9.0),
      (11.5, 10.0),
      (10.0, 11.5),
      (8.5, 10.0),
      (10.0, 8.5),
      (12.0, 10.0),
      (10.0, 12.0),
      (8.0, 10.0),
      (10.0, 8.0)
    ]
    for (p, t) in zip(ps, ts)
      @test p ≈ t
    end

    b = Ball(P3(0, 0, 0), T(2))
    ps = sample(b, RegularSampling(3, 2, 3))
    ts = P3[
      (0.7071067811865475, 0.0, 0.7071067811865476),
      (1.0, 0.0, 6.123233995736766e-17),
      (0.7071067811865476, 0.0, -0.7071067811865475),
      (-0.7071067811865475, 8.659560562354932e-17, 0.7071067811865476),
      (-1.0, 1.2246467991473532e-16, 6.123233995736766e-17),
      (-0.7071067811865476, 8.659560562354934e-17, -0.7071067811865475),
      (1.0606601717798212, 0.0, 1.0606601717798214),
      (1.5, 0.0, 9.184850993605148e-17),
      (1.0606601717798214, 0.0, -1.0606601717798212),
      (-1.0606601717798212, 1.2989340843532398e-16, 1.0606601717798214),
      (-1.5, 1.8369701987210297e-16, 9.184850993605148e-17),
      (-1.0606601717798214, 1.29893408435324e-16, -1.0606601717798212),
      (1.414213562373095, 0.0, 1.4142135623730951),
      (2.0, 0.0, 1.2246467991473532e-16),
      (1.4142135623730951, 0.0, -1.414213562373095),
      (-1.414213562373095, 1.7319121124709863e-16, 1.4142135623730951),
      (-2.0, 2.4492935982947064e-16, 1.2246467991473532e-16),
      (-1.4142135623730951, 1.7319121124709868e-16, -1.414213562373095)
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
    @test length(cs) == 200
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
    @test length(cs) == 200
    @test all(T(0) - eps(T) ≤ z ≤ T(1) + eps(T) for z in zs)

    # cylinder surface with non-parallel planes
    c = CylinderSurface(Plane(P3(0, 0, 0), V3(1, 0, 1)), Plane(P3(1, 1, 1), V3(0, 1, 1)), T(1))
    ps = sample(c, RegularSampling(20, 10))
    cs = coordinates.(ps)
    @test length(cs) == 200

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

    h = Hexahedron(P3[(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)])
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

    grid = CartesianGrid{T}(10, 10)
    points = sample(grid, RegularSampling(100, 200))
    @test length(collect(points)) == 20000

    torus = Torus(P3(0, 0, 0), V3(1, 0, 0), T(2), T(1))
    ps = sample(torus, RegularSampling(3, 3))
    ts = P3[
      (0, 0, 1),
      (0, -0.8660254037844387, -0.5),
      (0, 0.8660254037844387, -0.5),
      (-1, 0, 2),
      (-1, -1.7320508075688774, -1),
      (-1, 1.7320508075688774, -1),
      (1, 0, 2),
      (1, -1.7320508075688774, -1),
      (1, 1.7320508075688774, -1)
    ]
    for (p, t) in zip(ps, ts)
      @test p ≈ t
    end
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
    poly = PolyArea(
      P2[
        (-44.20065308, -21.12284851),
        (-44.20324135, -21.122799875),
        (-44.20582962, -21.12275124),
      ]
    )
    ps = sample(poly, MinDistanceSampling(3.2423333333753135e-5))
    @test length(ps) > 0
  end

  @testset "Utilities" begin
    # uniform sampling
    d = CartesianGrid{T}(10, 10)
    s = sample(d, 50)
    @test nelements(s) == 50
    @test s[1] isa Quadrangle

    # weighted sampling
    d = CartesianGrid{T}(10, 10, 10)
    s = sample(d, 100, rand([1, 2], 1000))
    @test nelements(s) == 100
    @test s[1] isa Hexahedron

    # ordered sampling
    d = CartesianGrid{T}(10, 10, 10)
    s = sample(d, 100, rand([1, 2], 1000), ordered=true)
    @test nelements(s) == 100
    @test s[1] isa Hexahedron
  end

  @testset "RNGs" begin
    dom = CartesianGrid{T}(100, 100)
    for method in [UniformSampling(100), WeightedSampling(100), BallSampling(T(10))]
      rng = MersenneTwister(2021)
      s1 = sample(rng, dom, method)
      rng = MersenneTwister(2021)
      s2 = sample(rng, dom, method)
      @test collect(s1) == collect(s2)
    end

    # cannot test some sampling methods with T = Float32
    # because of https://github.com/JuliaStats/StatsBase.jl/issues/695
    if T == Float64
      for method in [HomogeneousSampling(100), MinDistanceSampling(T(5))]
        rng = MersenneTwister(2021)
        s1 = sample(rng, dom, method)
        rng = MersenneTwister(2021)
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
      Hexahedron(P3[(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)])
    ]
      rng = MersenneTwister(2021)
      s1 = sample(rng, geom, method)
      rng = MersenneTwister(2021)
      s2 = sample(rng, geom, method)
      @test collect(s1) == collect(s2)
    end
  end
end
