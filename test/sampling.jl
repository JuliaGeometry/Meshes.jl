@testset "Sampling" begin
  @testset "RegularSampling" begin
    b = Box(P2(0, 0), P2(2, 2))
    ps = sample(b, RegularSampling(3))
    @test collect(ps) == P2[(0,0),(1,0),(2,0),(0,1),(1,1),(2,1),(0,2),(1,2),(2,2)]
    ps = sample(b, RegularSampling(2, 3))
    @test collect(ps) == P2[(0,0),(2,0),(0,1),(2,1),(0,2),(2,2)]

    s = Sphere(P2(0, 0), T(2))
    ps = sample(s, RegularSampling(4))
    ts = P2[(2,0),(0,2),(-2,0),(0,-2)]
    for (p, t) in zip(ps, ts)
      @test p ≈ t
    end

    s = Sphere(P3(0, 0, 0), T(2))
    ps = sample(s, RegularSampling(2,2))
    ts = P3[(1.7320508075688772, 0.0, 1.0),
            (1.7320508075688772, 0.0, -1.0),
            (-1.7320508075688772, 0.0, 1.0),
            (-1.7320508075688772, 0.0, -1.0)]
    for (p, t) in zip(ps, ts)
      @test p ≈ t
    end

    b = Ball(P2(0, 0), T(2))
    ps = sample(b, RegularSampling(4,3))
    ts = P2[(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0),
            (1.5, 0.0), (0.0, 1.5), (-1.5, 0.0), (0.0, -1.5),
            (2.0, 0.0), (0.0, 2.0), (-2.0, 0.0), (0.0, -2.0)]
    for (p, t) in zip(ps, ts)
      @test p ≈ t
    end

    b = Ball(P3(0, 0, 0), T(2))
    ps = sample(b, RegularSampling(3,2,3))
    ts = P3[(0.7071067811865475, 0.0, 0.7071067811865476),
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
            (-1.4142135623730951, 1.7319121124709868e-16, -1.414213562373095)]
    for (p, t) in zip(ps, ts)
      @test p ≈ t
    end
  end

  @testset "UniformSampling" begin
    Random.seed!(2021)
    d = CartesianGrid{T}(100,100)
    s = sample(d, UniformSampling(100))
    μ = mean(coordinates(s, 1:nelements(s)), dims=2)
    @test nelements(s) == 100
    @test isapprox(μ, [50.,50.], atol=T(10))
  end

  @testset "WeightedSampling" begin
    # uniform weights => uniform sampler
    Random.seed!(2020)
    d = CartesianGrid{T}(100,100)
    s = sample(d, WeightedSampling(100))
    μ = mean(coordinates(s, 1:nelements(s)), dims=2)
    @test nelements(s) == 100
    @test isapprox(μ, [50.,50.], atol=T(10))
  end

  @testset "BallSampling" begin
    d = CartesianGrid{T}(100,100)
    s = sample(d, BallSampling(T(10)))
    n = nelements(s)
    X = coordinates(s, sample(1:n, 2, replace=false))
    x, y = X[:,1], X[:,2]
    @test n < 100
    @test sqrt(sum((x - y).^2)) ≥ T(10)

    d = CartesianGrid{T}(100,100)
    s = sample(d, BallSampling(T(20)))
    n = nelements(s)
    X = coordinates(s, sample(1:n, 2, replace=false))
    x, y = X[:,1], X[:,2]
    @test n < 50
    @test sqrt(sum((x - y).^2)) ≥ T(20)
  end

  @testset "Utility" begin
    # uniform sampling
    d = CartesianGrid{T}(10,10)
    s = sample(d, 50)
    @test nelements(s) == 50
    @test s[1] isa Quadrangle

    # weighted sampling
    d = CartesianGrid{T}(10,10,10)
    s = sample(d, 100, rand([1,2], 1000))
    @test nelements(s) == 100
    @test s[1] isa Hexahedron
  end
end
