@testset "Neighborhoods" begin
  @testset "NormBall" begin
    # Euclidean metric
    ball = NormBall(T(1/2))
    @test evaluate(metric(ball), T[0], T[0]) ≤ radius(ball)
    @test evaluate(metric(ball), T[0], T[1]) > radius(ball)

    ball = NormBall(T(1))
    @test evaluate(metric(ball), T[0,0], T[0,0]) ≤ radius(ball)
    @test evaluate(metric(ball), T[0,0], T[1,0]) ≤ radius(ball)
    @test evaluate(metric(ball), T[0,0], T[0,1]) ≤ radius(ball)

    # Chebyshev metric
    ball = NormBall(T(1/2), Chebyshev())
    @test evaluate(metric(ball), T[0], T[0]) ≤ radius(ball)
    @test evaluate(metric(ball), T[0], T[1]) > radius(ball)

    for r in [1.,2.,3.,4.,5.]
      ball = NormBall(r, Chebyshev())
      for i in 0.0:1.0:r, j in 0.0:1.0:r
        @test evaluate(metric(ball), T[0,0], T[i,j]) ≤ radius(ball)
      end
    end

    ball = NormBall(T(1))
    if T == Float32
      @test sprint(show, ball) == "NormBall{Float32}(1.0, Euclidean(0.0))"
    else
      @test sprint(show, ball) == "NormBall{Float64}(1.0, Euclidean(0.0))"
    end
  end

  @testset "Ellipsoid" begin
    # 2D ellipsoid rotated 45 degrees in GSLIB convention
    ellipsoid = Ellipsoid(T[2,1], T[45], convention=GSLIB)

    # tests along main semiaxes, slightly below threshold
    @test evaluate(metric(ellipsoid), T[0,0], T(1.9) * T[√2/2, √2/2]) ≤ T(1)
    @test evaluate(metric(ellipsoid), T[0,0], T(0.9) * T[√2/2,-√2/2]) ≤ T(1)

    # tests along main semiaxes, slightly above threshold
    @test evaluate(metric(ellipsoid), T[0,0], T(2.1) * T[√2/2, √2/2]) > T(1)
    @test evaluate(metric(ellipsoid), T[0,0], T(1.1) * T[√2/2,-√2/2]) > T(1)

    # 3D ellipsoid rotated (45, -45, 0) in GSLIB convention
    ellipsoid = Ellipsoid(T[3,2,1], T[45,-45,0], convention=GSLIB)

    # tests along main semiaxes, slightly below threshold
    @test evaluate(metric(ellipsoid), T[0,0,0], T(2.9) * T[0.5,0.5,-√2/2]) ≤ T(1)
    @test evaluate(metric(ellipsoid), T[0,0,0], T(1.9) * T[√2/2,-√2/2,0.0]) ≤ T(1)
    @test evaluate(metric(ellipsoid), T[0,0,0], T(0.9) * T[0.5,0.5,√2/2]) ≤ T(1)

    # Tests along main semiaxes, slightly above threshold
    @test evaluate(metric(ellipsoid), T[0,0,0], T(3.1) * T[0.5,0.5,-√2/2]) > T(1)
    @test evaluate(metric(ellipsoid), T[0,0,0], T(2.1) * T[√2/2,-√2/2,0.0]) > T(1)
    @test evaluate(metric(ellipsoid), T[0,0,0], T(1.1) * T[0.5,0.5,√2/2]) > T(1)

    ellipsoid = Ellipsoid(T.((2,1)), T.((45,)), convention=GSLIB)
    if T == Float32
      @test sprint(show, ellipsoid) == "Ellipsoid{Float32}((2.0f0, 1.0f0), (45.0f0,), GSLIB)"
    else
      @test sprint(show, ellipsoid) == "Ellipsoid{Float64}((2.0, 1.0), (45.0,), GSLIB)"
    end
  end
end
