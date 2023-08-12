@testset "Neighbor search" begin
  @testset "BallSearch" begin
    𝒟 = CartesianGrid((10, 10), T.((-0.5, -0.5)), T.((1.0, 1.0)))

    S = BallSearch(𝒟, MetricBall(T(1)))
    n = search(P2(0, 0), S)
    @test Set(n) == Set([1, 2, 11])
    n = search(P2(9, 0), S)
    @test Set(n) == Set([9, 10, 20])
    n = search(P2(0, 9), S)
    @test Set(n) == Set([91, 81, 92])
    n = search(P2(9, 9), S)
    @test Set(n) == Set([100, 99, 90])

    S = BallSearch(𝒟, MetricBall(T(√2 + eps(T))))
    n = search(P2(0, 0), S)
    @test Set(n) == Set([1, 2, 11, 12])
    n = search(P2(9, 0), S)
    @test Set(n) == Set([9, 10, 19, 20])
    n = search(P2(0, 9), S)
    @test Set(n) == Set([81, 82, 91, 92])
    n = search(P2(9, 9), S)
    @test Set(n) == Set([89, 90, 99, 100])

    # non MinkowskiMetric example
    𝒟 = CartesianGrid((360, 180), T.((0.0, -90.0)), T.((1.0, 1.0)))
    S = BallSearch(𝒟, MetricBall(T(150), Haversine(T(6371))))
    n = search(P2(0, 0), S)
    @test Set(n) == Set([32041, 32400, 32401, 32760])
  end

  @testset "KNearestSearch" begin
    𝒟 = CartesianGrid((10, 10), T.((-0.5, -0.5)), T.((1.0, 1.0)))
    s = KNearestSearch(𝒟, 3)
    n = search(P2(0, 0), s)
    @test Set(n) == Set([1, 2, 11])
    n = search(P2(9, 0), s)
    @test Set(n) == Set([9, 10, 20])
    n = search(P2(0, 9), s)
    @test Set(n) == Set([91, 81, 92])
    n = search(P2(9, 9), s)
    @test Set(n) == Set([100, 99, 90])
    n, d = searchdists(P2(9, 9), s)
    @test Set(n) == Set([100, 99, 90])
    @test length(d) == 3
    n = Vector{Int}(undef, maxneighbors(s))
    nn = search!(n, P2(9, 9), s)
    @test nn == 3
    @test Set(n[1:nn]) == Set([100, 99, 90])
    n = Vector{Int}(undef, maxneighbors(s))
    d = Vector{T}(undef, maxneighbors(s))
    nn = searchdists!(n, d, P2(9, 9), s)
    @test nn == 3
    @test Set(n[1:nn]) == Set([100, 99, 90])

    𝒟 = CartesianGrid((10, 10), T.((-0.5, -0.5)), T.((1.0, 1.0)))
    s = KNearestSearch(𝒟, 3)
    skipfun(i) = i == 1 || i == 2 || i == 11
    n = search(P2(0, 0), s, skip=skipfun)
    @test Set(n) == Set([12, 21, 3])
  end

  @testset "KBallSearch" begin
    𝒟 = CartesianGrid((10, 10), T.((-0.5, -0.5)), T.((1.0, 1.0)))

    s = KBallSearch(𝒟, 10, MetricBall(T(100)))
    n = search(P2(5, 5), s)
    @test length(n) == 10

    s = KBallSearch(𝒟, 10, MetricBall(T.((100, 100))))
    n = search(P2(5, 5), s)
    @test length(n) == 10

    s = KBallSearch(𝒟, 10, MetricBall(T(1)))
    n = search(P2(5, 5), s)
    @test length(n) == 5
    @test n[1] == 56

    s = KBallSearch(𝒟, 10, MetricBall(T(1)))
    n, d = searchdists(P2(5, 5), s)
    @test length(n) == 5
    @test length(d) == 5

    s = KBallSearch(𝒟, 10, MetricBall(T(1)))
    n = Vector{Int}(undef, maxneighbors(s))
    nn = search!(n, P2(5, 5), s)
    @test nn == 5

    s = KBallSearch(𝒟, 10, MetricBall(T(1)))
    n = Vector{Int}(undef, maxneighbors(s))
    d = Vector{T}(undef, maxneighbors(s))
    nn = searchdists!(n, d, P2(5, 5), s)
    @test nn == 5

    s = KBallSearch(𝒟, 10, MetricBall(T(1)))
    skipfun(i) = i == 56
    n = search(P2(5, 5), s, skip=skipfun)
    @test length(n) == 4
    n = search(P2(-0.2, -0.2), s)
    @test length(n) == 1
    n = search(P2(-10, -10), s)
    @test length(n) == 0
    n, d = searchdists(P2(5, 5), s, skip=skipfun)
    @test length(n) == 4
    @test length(d) == 4
  end
end
