@testitem "MetricBall" setup = [Setup] begin
  # Euclidean metric
  b = MetricBall(T(1 / 2))
  r = radius(b)
  m = metric(b)
  @test m(cart(0, 0), cart(0, 0)) ≤ r
  @test m(cart(0, 0), cart(1, 0)) > r
  @test radii(b) == (T(1 / 2) * u"m",)
  b = MetricBall(T(1))
  r = radius(b)
  m = metric(b)
  @test m(cart(0, 0), cart(0, 0)) ≤ r
  @test m(cart(0, 0), cart(1, 0)) ≤ r
  @test m(cart(0, 0), cart(0, 1)) ≤ r
  @test hasequalradii(b)
  if T === Float32
    @test sprint(show, b) == "MetricBall(1.0f0 m, EuclideanDistance)"
  else
    @test sprint(show, b) == "MetricBall(1.0 m, EuclideanDistance)"
  end

  # test radius convention
  b = MetricBall(T.((1, 1)))
  m = metric(b)
  @test radius(b) == oneunit(ℳ)
  @test m(cart(1, 0), cart(0, 0)) == m(cart(0, 1), cart(0, 0))
  b = MetricBall(T.((1, 2)))
  m = metric(b)
  @test radius(b) == oneunit(ℳ)
  @test m(cart(1, 0), cart(0, 0)) != m(cart(0, 1), cart(0, 0))
  b = MetricBall(T.((1.0, 0.5, 0.5)), RotZYX(T(-π / 4), T(0), T(0)))
  m = metric(b)
  @test radius(b) == oneunit(ℳ)
  @test m(cart(1, 1, 0), cart(0, 0, 0)) ≈ √T(8) * u"m"
  @test m(cart(-1, 1, 0), cart(0, 0, 0)) ≈ √T(2) * u"m"

  # make sure the correct constructor is called
  m = metric(MetricBall(T.((1.0, 0.5, 0.2)), RotXYX(T(0), T(0), T(0))))
  @test m isa MahalanobisDistance

  # make sure the angle is clockwise
  m = metric(MetricBall(T.((20.0, 5.0)), Angle2d(T(π / 2))))
  @test m isa MahalanobisDistance
  @test m(cart(1, 0), cart(0, 0)) ≈ T(0.2) * u"m"
  @test m(cart(0, 1), cart(0, 0)) ≈ T(0.05) * u"m"

  # basic multiplication
  @test 2MetricBall(T(1)) == MetricBall(T(2))
  @test 2MetricBall(T.((1, 2, 3))) == MetricBall(T.((2, 4, 6)))

  # access to rotation
  @test rotation(MetricBall(T(1))) == I
  @test rotation(MetricBall(T.((1, 2, 3)))) == I
  @test rotation(MetricBall(T.((1, 2)), Angle2d(T(π / 2)))) == Angle2d(T(π / 2))

  # make sure radii are converted to unitful float
  @test MetricBall(1) == MetricBall(1.0u"m")
  @test MetricBall(1) == MetricBall(1.0u"m")
end
