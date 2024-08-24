@testitem "MetricBall" setup = [Setup] begin
  # Euclidean metric
  b = MetricBall(T(1 / 2))
  r = radius(b)
  m = metric(b)
  @test evaluate(m, T[0] * u"m", T[0] * u"m") ≤ r
  @test evaluate(m, T[0] * u"m", T[1] * u"m") > r
  @test radii(b) == (T(1 / 2) * u"m",)

  b = MetricBall(T(1))
  r = radius(b)
  m = metric(b)
  @test evaluate(m, T[0, 0] * u"m", T[0, 0] * u"m") ≤ r
  @test evaluate(m, T[0, 0] * u"m", T[1, 0] * u"m") ≤ r
  @test evaluate(m, T[0, 0] * u"m", T[0, 1] * u"m") ≤ r
  @test isisotropic(b)
  if T === Float32
    @test sprint(show, b) == "MetricBall(1.0f0 m, Euclidean)"
  else
    @test sprint(show, b) == "MetricBall(1.0 m, Euclidean)"
  end

  # Chebyshev metric
  b = MetricBall(T(1 / 2), Chebyshev())
  r = radius(b)
  m = metric(b)
  @test evaluate(m, T[0] * u"m", T[0] * u"m") ≤ r
  @test evaluate(m, T[0] * u"m", T[1] * u"m") > r

  for r in T[1.0, 2.0, 3.0, 4.0, 5.0]
    b = MetricBall(r, Chebyshev())
    r = radius(b)
    m = metric(b)
    for i in zero(r):oneunit(r):r, j in zero(r):oneunit(r):r
      @test evaluate(m, T[0, 0] * u"m", [i, j]) ≤ r
    end
  end

  # 2D simple test of default convention
  b = MetricBall(T.((1, 1)))
  m = metric(b)
  @test radius(b) == oneunit(ℳ)
  @test evaluate(m, T[1, 0] * u"m", T[0, 0] * u"m") == evaluate(m, T[0, 1] * u"m", T[0, 0] * u"m")

  b = MetricBall(T.((1, 2)))
  m = metric(b)
  @test radius(b) == oneunit(ℳ)
  @test evaluate(m, T[1, 0] * u"m", T[0, 0] * u"m") != evaluate(m, T[0, 1] * u"m", T[0, 0] * u"m")

  # 3D simple test of default convention
  b = MetricBall(T.((1.0, 0.5, 0.5)), RotZYX(T(-π / 4), T(0), T(0)))
  m = metric(b)
  @test radius(b) == oneunit(ℳ)
  @test evaluate(m, T[1.0, 1.0, 0.0] * u"m", T[0.0, 0.0, 0.0] * u"m") ≈ √T(8) * u"m"
  @test evaluate(m, T[-1.0, 1.0, 0.0] * u"m", T[0.0, 0.0, 0.0] * u"m") ≈ √T(2) * u"m"

  # make sure the correct constructor is called
  m = metric(MetricBall(T.((1.0, 0.5, 0.2)), RotXYX(T(0), T(0), T(0))))
  @test m isa Mahalanobis

  # make sure the angle is clockwise
  m = metric(MetricBall(T.((20.0, 5.0)), Angle2d(T(π / 2))))
  @test m isa Mahalanobis
  @test evaluate(m, T[1.0, 0.0] * u"m", T[0.0, 0.0] * u"m") ≈ T(0.2) * u"m"
  @test evaluate(m, T[0.0, 1.0] * u"m", T[0.0, 0.0] * u"m") ≈ T(0.05) * u"m"

  # basic multiplication
  @test 2MetricBall(T(1)) == MetricBall(T(2))
  @test 2MetricBall(T.((1, 2, 3))) == MetricBall(T.((2, 4, 6)))

  # access to rotation
  @test rotation(MetricBall(T(1))) == I
  @test rotation(MetricBall(T.((1, 2, 3)))) == I
  @test rotation(MetricBall(T.((1, 2)), Angle2d(T(π / 2)))) == Angle2d(T(π / 2))
end
