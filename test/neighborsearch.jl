@testitem "BallSearch" setup = [Setup] begin
  𝒟 = CartesianGrid(T.((-0.5, -0.5)), T.((1.0, 1.0)), GridTopology(10, 10))

  s = BallSearch(𝒟, MetricBall(T(1)))
  n = search(cart(0, 0), s)
  @test Set(n) == Set([1, 2, 11])
  n = search(cart(9, 0), s)
  @test Set(n) == Set([9, 10, 20])
  n = search(cart(0, 9), s)
  @test Set(n) == Set([91, 81, 92])
  n = search(cart(9, 9), s)
  @test Set(n) == Set([100, 99, 90])

  s = BallSearch(𝒟, MetricBall(T(√2 + eps(T))))
  n = search(cart(0, 0), s)
  @test Set(n) == Set([1, 2, 11, 12])
  n = search(cart(9, 0), s)
  @test Set(n) == Set([9, 10, 19, 20])
  n = search(cart(0, 9), s)
  @test Set(n) == Set([81, 82, 91, 92])
  n = search(cart(9, 9), s)
  @test Set(n) == Set([89, 90, 99, 100])

  # different units
  s = BallSearch(𝒟, MetricBall(T(10) * u"dm"))
  n = search(Point(T(900) * u"cm", T(900) * u"cm"), s)
  @test Set(n) == Set([100, 99, 90])
  n = search(Point(T(9000) * u"mm", T(9000) * u"mm"), s)
  @test Set(n) == Set([100, 99, 90])

  # construct from vector of geometries
  s = BallSearch([cart(1, 1), cart(2, 2), cart(3, 3)], MetricBall(T(1)))
  @test s isa BallSearch
end

@testitem "KNearestSearch" setup = [Setup] begin
  𝒟 = CartesianGrid(T.((-0.5, -0.5)), T.((1.0, 1.0)), GridTopology(10, 10))
  s = KNearestSearch(𝒟, 3)
  n = search(cart(0, 0), s)
  @test Set(n) == Set([1, 2, 11])
  n = search(cart(9, 0), s)
  @test Set(n) == Set([9, 10, 20])
  n = search(cart(0, 9), s)
  @test Set(n) == Set([91, 81, 92])
  n = search(cart(9, 9), s)
  @test Set(n) == Set([100, 99, 90])
  n, d = searchdists(cart(9, 9), s)
  @test Set(n) == Set([100, 99, 90])
  @test length(d) == 3
  n = Vector{Int}(undef, maxneighbors(s))
  nn = search!(n, cart(9, 9), s)
  @test nn == 3
  @test Set(n[1:nn]) == Set([100, 99, 90])
  n = Vector{Int}(undef, maxneighbors(s))
  d = Vector{ℳ}(undef, maxneighbors(s))
  nn = searchdists!(n, d, cart(9, 9), s)
  @test nn == 3
  @test Set(n[1:nn]) == Set([100, 99, 90])

  # different units
  s = KNearestSearch(𝒟, 3)
  n = search(Point(T(900) * u"cm", T(900) * u"cm"), s)
  @test Set(n) == Set([100, 99, 90])
  n = search(Point(T(9000) * u"mm", T(9000) * u"mm"), s)
  @test Set(n) == Set([100, 99, 90])

  # construct from vector of geometries
  s = KNearestSearch([cart(1, 1), cart(2, 2), cart(3, 3)], 3)
  @test s isa KNearestSearch
end

@testitem "KBallSearch" setup = [Setup] begin
  𝒟 = CartesianGrid(T.((-0.5, -0.5)), T.((1.0, 1.0)), GridTopology(10, 10))

  s = KBallSearch(𝒟, 10, MetricBall(T(100)))
  n = search(cart(5, 5), s)
  @test length(n) == 10

  s = KBallSearch(𝒟, 10, MetricBall(T.((100, 100))))
  n = search(cart(5, 5), s)
  @test length(n) == 10

  s = KBallSearch(𝒟, 10, MetricBall(T(1)))
  n = search(cart(5, 5), s)
  @test length(n) == 5
  @test n[1] == 56

  s = KBallSearch(𝒟, 10, MetricBall(T(1)))
  n, d = searchdists(cart(5, 5), s)
  @test length(n) == 5
  @test length(d) == 5

  s = KBallSearch(𝒟, 10, MetricBall(T(1)))
  n = Vector{Int}(undef, maxneighbors(s))
  nn = search!(n, cart(5, 5), s)
  @test nn == 5

  s = KBallSearch(𝒟, 10, MetricBall(T(1)))
  n = Vector{Int}(undef, maxneighbors(s))
  d = Vector{ℳ}(undef, maxneighbors(s))
  nn = searchdists!(n, d, cart(5, 5), s)
  @test nn == 5

  mask = trues(nelements(𝒟))
  mask[56] = false
  n = search(cart(5, 5), s, mask=mask)
  @test length(n) == 4
  n = search(cart(-0.2, -0.2), s)
  @test length(n) == 1
  n = search(cart(-10, -10), s)
  @test length(n) == 0
  n, d = searchdists(cart(5, 5), s, mask=mask)
  @test length(n) == 4
  @test length(d) == 4

  # different units
  s = KBallSearch(𝒟, 10, MetricBall(T(10) * u"dm"))
  n = search(Point(T(500) * u"cm", T(500) * u"cm"), s)
  @test Set(n) == Set([56, 66, 55, 57, 46])
  n = search(Point(T(5000) * u"mm", T(5000) * u"mm"), s)
  @test Set(n) == Set([56, 66, 55, 57, 46])

  # construct from vector of geometries
  s = KBallSearch([cart(1, 1), cart(2, 2), cart(3, 3)], 3, MetricBall(T(1)))
  @test s isa KBallSearch
end
