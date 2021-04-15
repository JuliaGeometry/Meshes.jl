@testset "Neighborhood search" begin
  @testset "NeighborhoodSearch" begin
    𝒟 = CartesianGrid((10, 10), T.((-0.5,-0.5)), T.((1.0,1.0)))

    S = NeighborhoodSearch(𝒟, NormBall(T(1)))
    n = search(P2(0,0), S)
    @test Set(n) == Set([1,2,11])
    n = search(P2(9,0), S)
    @test Set(n) == Set([9,10,20])
    n = search(P2(0,9), S)
    @test Set(n) == Set([91,81,92])
    n = search(P2(9,9), S)
    @test Set(n) == Set([100,99,90])

    S = NeighborhoodSearch(𝒟, NormBall(T(√2+eps(T))))
    n = search(P2(0,0), S)
    @test Set(n) == Set([1,2,11,12])
    n = search(P2(9,0), S)
    @test Set(n) == Set([9,10,19,20])
    n = search(P2(0,9), S)
    @test Set(n) == Set([81,82,91,92])
    n = search(P2(9,9), S)
    @test Set(n) == Set([89,90,99,100])

    # non MinkowskiMetric example
    𝒟 = CartesianGrid((360,180), T.((0.0,-90.0)), T.((1.0,1.0)))
    S = NeighborhoodSearch(𝒟, NormBall(T(150), Haversine(T(6371))))
    n = search(P2(0,0), S)
    @test Set(n) == Set([32041, 32400, 32401, 32760])
  end

  @testset "KNearestSearch" begin
    𝒟 = CartesianGrid((10,10), T.((-0.5,-0.5)), T.((1.0,1.0)))
    S = KNearestSearch(𝒟, 3)
    n = search(P2(0,0), S)
    @test Set(n) == Set([1,2,11])
    n = search(P2(9,0), S)
    @test Set(n) == Set([9,10,20])
    n = search(P2(0,9), S)
    @test Set(n) == Set([91,81,92])
    n = search(P2(9,9), S)
    @test Set(n) == Set([100,99,90])
  end

  @testset "KBallSearch" begin
    𝒟 = CartesianGrid((10,10), T.((-0.5,-0.5)), T.((1.0,1.0)))

    s = KBallSearch(𝒟, 10, NormBall(T(100)))
    n = search(P2(5,5), s)
    @test length(n) == 10

    s = KBallSearch(𝒟, 10, NormBall(T(1)))
    n = search(P2(5,5), s)
    @test length(n) == 5
    @test n[1] == 56

    mask = trues(nelements(𝒟))
    mask[56] = false
    n = search(P2(5,5), s, mask=mask)
    @test length(n) == 4
    n = search(P2(-0.2,-0.2), s)
    @test length(n) == 1
    n = search(P2(-10,-10), s)
    @test length(n) == 0
  end

  @testset "FilteredSearch" begin
	# dummy type implementing the Data trait
	struct DummyData{𝒟,𝒯} <: Data
		domain::𝒟
		table::𝒯
	end
	Meshes.domain(data::DummyData) = data.domain
	Meshes.values(data::DummyData) = data.table

	# nmax tests
	𝒟 = CartesianGrid((10,10), T.((-0.5,-0.5)), T.((1.0,1.0)))
	S1 = NeighborhoodSearch(𝒟, NormBall(T(5)))
	S2 = KNearestSearch(𝒟, 10)
	B1 = FilteredSearch(S1, 5)
	B2 = FilteredSearch(S2, 5)
	p = centroid(𝒟, rand(1:100))
	n = search(p, B1)
	@test length(n) == 5
	p = centroid(𝒟, rand(1:100))
	n = search(p, B2)
	@test length(n) == 5

	# maxpercategory tests
	catgs = rand(1:2,100)
	props = (a=catgs, b=catgs)
	data  = DummyData(𝒟, props)
	S3 = NeighborhoodSearch(data, Ellipsoid([4,2],[pi/4]))
	B3 = FilteredSearch(S3, maxpercategory=(a=5,b=5))
	p  = Point([4.5, 4.5])
	n  = search(p, B3)
	c  = view(catgs, n)
	@test sum(c .== 1) == 5
	@test sum(c .== 2) == 5
	@test length(n) == 10

	# maxpersector tests
	𝒟  = CartesianGrid((10,10,10), T.((-5,-5,-5)), T.((1.0,1.0,1.0)))
	S4 = NeighborhoodSearch(𝒟, Ellipsoid([10,10,10],[0.001,0.001,0.001]))
	B4 = FilteredSearch(S4, 30, maxpersector=3)
	p = Point([0, 0, 0])
	n = search(p, B4)
	@test length(n) == 24
	coords = coordinates.(centroid.(view(𝒟, n)))
	@test sum([all(>(0), x) for x in coords]) == 3
	@test sum([all(<(0), x) for x in coords]) == 3
  end
end
