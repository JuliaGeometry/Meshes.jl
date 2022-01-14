@testset "Neighborhoods" begin
  @testset "MetricBall" begin
    # Euclidean metric
    b = MetricBall(T(1/2))
    r = radius(b)
    m = metric(b)
    @test evaluate(m, T[0], T[0]) ≤ r
    @test evaluate(m, T[0], T[1]) > r

    b = MetricBall(T(1))
    r = radius(b)
    m = metric(b)
    @test evaluate(m, T[0,0], T[0,0]) ≤ r
    @test evaluate(m, T[0,0], T[1,0]) ≤ r
    @test evaluate(m, T[0,0], T[0,1]) ≤ r

    # Chebyshev metric
    b = MetricBall(T(1/2), Chebyshev())
    r = radius(b)
    m = metric(b)
    @test evaluate(m, T[0], T[0]) ≤ r
    @test evaluate(m, T[0], T[1]) > r

    for r in [1.,2.,3.,4.,5.]
      b = MetricBall(r, Chebyshev())
      r = radius(b)
      m = metric(b)
      for i in 0.0:1.0:r, j in 0.0:1.0:r
        @test evaluate(m, T[0,0], T[i,j]) ≤ r
      end
    end

    # 2D simple test of default convention
    m = metric(MetricBall(T.((1,1))))
    @test evaluate(m, T[1,0], T[0,0]) == evaluate(m, T[0,1], T[0,0])

    m = metric(MetricBall(T.((1,2))))
    @test evaluate(m, T[1,0], T[0,0]) != evaluate(m, T[0,1], T[0,0])

    # 3D simple test of default convention
    m = metric(MetricBall(T.((1.0,0.5,0.5)), TaitBryanAngles(T(π/4),T(0),T(0))))
    @test evaluate(m, [1.,1.,0.], [0.,0.,0.]) ≈ √T(2)
    @test evaluate(m, [-1.,1.,0.], [0.,0.,0.]) ≈ √T(8)
  end
end
