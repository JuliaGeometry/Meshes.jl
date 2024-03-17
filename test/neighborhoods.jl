@testset "Neighborhoods" begin
  @testset "MetricBall" begin
    # Euclidean metric
    b = MetricBall(T(1 / 2))
    r = radius(b)
    m = metric(b)
    @test evaluate(m, T[0], T[0]) ≤ r
    @test evaluate(m, T[0], T[1]) > r
    @test radii(b) == T[1 / 2]

    b = MetricBall(T(1))
    r = radius(b)
    m = metric(b)
    @test evaluate(m, T[0, 0], T[0, 0]) ≤ r
    @test evaluate(m, T[0, 0], T[1, 0]) ≤ r
    @test evaluate(m, T[0, 0], T[0, 1]) ≤ r
    @test isisotropic(b)
    @test sprint(show, b) == "MetricBall(1.0, Euclidean)"

    # Chebyshev metric
    b = MetricBall(T(1 / 2), Chebyshev())
    r = radius(b)
    m = metric(b)
    @test evaluate(m, T[0], T[0]) ≤ r
    @test evaluate(m, T[0], T[1]) > r

    for r in [1.0, 2.0, 3.0, 4.0, 5.0]
      b = MetricBall(r, Chebyshev())
      r = radius(b)
      m = metric(b)
      for i in 0.0:1.0:r, j in 0.0:1.0:r
        @test evaluate(m, T[0, 0], T[i, j]) ≤ r
      end
    end

    # 2D simple test of default convention
    m = metric(MetricBall(T.((1, 1))))
    @test evaluate(m, T[1, 0], T[0, 0]) == evaluate(m, T[0, 1], T[0, 0])

    m = metric(MetricBall(T.((1, 2))))
    @test evaluate(m, T[1, 0], T[0, 0]) != evaluate(m, T[0, 1], T[0, 0])

    # 3D simple test of default convention
    m = metric(MetricBall(T.((1.0, 0.5, 0.5)), RotZYX(T(-π / 4), T(0), T(0))))
    @test evaluate(m, [1.0, 1.0, 0.0], [0.0, 0.0, 0.0]) ≈ √T(8)
    @test evaluate(m, [-1.0, 1.0, 0.0], [0.0, 0.0, 0.0]) ≈ √T(2)

    # make sure the correct constructor is called
    m = metric(MetricBall(T[1.0, 0.5, 0.2], RotXYX(T(0), T(0), T(0))))
    @test m isa Mahalanobis

    # make sure the angle is clockwise
    m = metric(MetricBall(T[20.0, 5.0], Angle2d(T(π / 2))))
    @test m isa Mahalanobis
    @test evaluate(m, [1.0, 0.0], [0.0, 0.0]) ≈ T(0.2)
    @test evaluate(m, [0.0, 1.0], [0.0, 0.0]) ≈ T(0.05)

    # basic multiplication
    @test 2MetricBall(T(1)) == MetricBall(T(2))
    @test 2MetricBall(T[1, 2, 3]) == MetricBall(T[2, 4, 6])

    # access to rotation
    @test rotation(MetricBall(T(1))) == I
    @test rotation(MetricBall(T[1, 2, 3])) == I
    @test rotation(MetricBall(T[1, 2], Angle2d(T(π / 2)))) == Angle2d(T(π / 2))
  end
end
