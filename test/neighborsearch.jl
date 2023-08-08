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
    S = KNearestSearch(𝒟, 3)
    n = search(P2(0, 0), S)
    @test Set(n) == Set([1, 2, 11])
    n = search(P2(9, 0), S)
    @test Set(n) == Set([9, 10, 20])
    n = search(P2(0, 9), S)
    @test Set(n) == Set([91, 81, 92])
    n = search(P2(9, 9), S)
    @test Set(n) == Set([100, 99, 90])
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

    mask = trues(nelements(𝒟))
    mask[56] = false
    n = search(P2(5, 5), s, mask=mask)
    @test length(n) == 4
    n = search(P2(-0.2, -0.2), s)
    @test length(n) == 1
    n = search(P2(-10, -10), s)
    @test length(n) == 0
  end

  @testset "GlobalSearch" begin
    𝒟 = CartesianGrid(10, 10)
    S = GlobalSearch(𝒟)
    p = centroid(𝒟, rand(1:100))
    n = search(p, S)
    @test n == 1:100
    mask = falses(nelements(𝒟))
    mask[15] = true
    mask[50] = true
    mask[90] = true
    n = search(p, S, mask=mask)
    @test n == [15, 50, 90]
  end
end
