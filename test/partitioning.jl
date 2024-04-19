@testset "Partitioning" begin
  setify(lists) = Set(Set.(lists))

  g = CartesianGrid{T}(10, 10)
  p = partition(g, UniformPartition(100))
  @test parent(p) == g
  @test length(p) == 100

  @testset "UniformPartition" begin
    rng = StableRNG(123)
    g = CartesianGrid{T}(3, 3)
    p = partition(rng, g, UniformPartition(3, false))
    @test setify(indices(p)) == setify([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    rng = StableRNG(123)
    p = partition(rng, g, UniformPartition(3))
    @test setify(indices(p)) == setify([[5, 4, 2], [6, 7, 8], [9, 3, 1]])

    g = CartesianGrid{T}(2, 3)
    p = partition(g, UniformPartition(3, false))
    @test setify(indices(p)) == setify([[1, 2], [3, 4], [5, 6]])

    # reproducible results with rng
    rng = StableRNG(123)
    g = CartesianGrid{T}(10, 10)
    p1 = partition(rng, g, UniformPartition(3))
    rng = StableRNG(123)
    p2 = partition(rng, g, UniformPartition(3))
    @test p1 == p2
  end

  @testset "DirectionPartition" begin
    g = CartesianGrid{T}(3, 3)

    # basic checks on small grids
    p = partition(g, DirectionPartition(T.((1, 0))))
    @test setify(indices(p)) == setify([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

    p = partition(g, DirectionPartition(T.((0, 1))))
    @test setify(indices(p)) == setify([[1, 4, 7], [2, 5, 8], [3, 6, 9]])

    p = partition(g, DirectionPartition(T.((1, 1))))
    @test setify(indices(p)) == setify([[1, 5, 9], [2, 6], [3], [4, 8], [7]])

    p = partition(g, DirectionPartition(T.((1, -1))))
    @test setify(indices(p)) == setify([[1], [2, 4], [3, 5, 7], [6, 8], [9]])

    # opposite directions produce same partition
    dir1 = (rand(T), rand(T))
    dir2 = .-dir1
    p1 = partition(g, DirectionPartition(dir1))
    p2 = partition(g, DirectionPartition(dir2))
    @test setify(indices(p1)) == setify(indices(p2))

    # partition of arbitrarily large grid always
    # returns the "lines" and "columns"
    for n in [10, 100, 200]
      g = CartesianGrid{T}(n, n)

      p = partition(g, DirectionPartition(T.((1, 0))))
      @test setify(indices(p)) == setify([collect(((i - 1) * n + 1):(i * n)) for i in 1:n])
      ns = [nelements(d) for d in p]
      @test all(ns .== n)

      p = partition(g, DirectionPartition(T.((0, 1))))
      @test setify(indices(p)) == setify([collect(i:n:(n * n)) for i in 1:n])
      ns = [nelements(d) for d in p]
      @test all(ns .== n)
    end

    # reproducible results with rng
    rng = StableRNG(123)
    g = CartesianGrid{T}(10, 10)
    p1 = partition(rng, g, DirectionPartition(T.((1, 0))))
    rng = StableRNG(123)
    p2 = partition(rng, g, DirectionPartition(T.((1, 0))))
    @test p1 == p2
  end

  @testset "FractionPartition" begin
    g = CartesianGrid{T}(10, 10)

    p = partition(g, FractionPartition(T(0.5)))
    @test nelements(p[1]) == nelements(p[2]) == 50
    @test length(p) == 2

    p = partition(g, FractionPartition(T(0.7)))
    @test nelements(p[1]) == 70
    @test nelements(p[2]) == 30

    p = partition(g, FractionPartition(T(0.3)))
    @test nelements(p[1]) == 30
    @test nelements(p[2]) == 70

    # reproducible results with rng
    rng = StableRNG(123)
    g = CartesianGrid{T}(10, 10)
    p1 = partition(rng, g, FractionPartition(T(0.5)))
    rng = StableRNG(123)
    p2 = partition(rng, g, FractionPartition(T(0.5)))
    @test p1 == p2
  end

  @testset "BlockPartition" begin
    g = CartesianGrid{T}(10, 10)

    p = partition(g, BlockPartition(T(5), T(5)))
    @test length(p) == 4
    @test all(nelements.(p) .== 25)

    p = partition(g, BlockPartition(T(5), T(2)))
    @test length(p) == 12
    @test Set(nelements.(p)) == Set([5, 10])

    g = CartesianGrid{T}(50, 50, 50)

    p = partition(g, BlockPartition(T(1.0), T(1.0), T(1.0), neighbors=false))
    @test length(p) == 125000
    @test Set(nelements.(p)) == Set(1)
    @test metadata(p) == Dict{Any,Any}()

    p = partition(g, BlockPartition(T(5.0), T(5.0), T(5.0), neighbors=true))
    @test length(p) == 1000
    @test Set(nelements.(p)) == Set(125)
    n = metadata(p)[:neighbors]
    @test length(n) == length(p)
    @test all(0 .< length.(n) .< 27)

    # reproducible results with rng
    rng = StableRNG(123)
    g = CartesianGrid{T}(10, 10)
    p1 = partition(rng, g, BlockPartition(T(5), T(2)))
    rng = StableRNG(123)
    p2 = partition(rng, g, BlockPartition(T(5), T(2)))
    @test p1 == p2

    m1 = BlockPartition((T(5), T(2)))
    m2 = BlockPartition(T(5), T(2))
    m3 = BlockPartition((T(5), T(2)), neighbors=false)
    m4 = BlockPartition(T(5), T(2), neighbors=false)
    @test m1 == m2 == m3 == m4
    m1 = BlockPartition(T(1))
    m2 = BlockPartition(T(1), neighbors=false)
    @test m1 == m2
  end

  @testset "BisectPointPartition" begin
    g = CartesianGrid((10, 10), T.((-0.5, -0.5)), T.((1.0, 1.0)))

    p = partition(g, BisectPointPartition(T.((0.0, 1.0)), T.((5.0, 5.1))))
    p1, p2 = p[1], p[2]
    @test nelements(p1) == 60
    @test nelements(p2) == 40

    # all points in p1 are below those in p2
    pts1 = [centroid(p1, i) for i in 1:nelements(p1)]
    pts2 = [centroid(p2, i) for i in 1:nelements(p2)]
    X1 = reduce(hcat, coordinates.(pts1))
    X2 = reduce(hcat, coordinates.(pts2))
    M1 = maximum(X1, dims=2)
    m2 = minimum(X2, dims=2)
    @test all(X1[2, j] < m2[2] for j in 1:size(X1, 2))
    @test all(X2[2, j] > M1[2] for j in 1:size(X2, 2))

    # flipping normal direction is equivalent to swapping subsets
    p₁ = partition(g, BisectPointPartition(T.((1.0, 0.0)), T.((5.1, 5.0))))
    p₂ = partition(g, BisectPointPartition(T.((-1.0, 0.0)), T.((5.1, 5.0))))
    @test nelements(p₁[1]) == nelements(p₂[2]) == 60
    @test nelements(p₁[2]) == nelements(p₂[1]) == 40

    # reproducible results with rng
    rng = StableRNG(123)
    g = CartesianGrid{T}(10, 10)
    p1 = partition(rng, g, BisectPointPartition(T.((1, 0)), T.((5, 5))))
    rng = StableRNG(123)
    p2 = partition(rng, g, BisectPointPartition(T.((1, 0)), T.((5, 5))))
    @test p1 == p2
  end

  @testset "BisectFractionPartition" begin
    g = CartesianGrid((10, 10), T.((-0.5, -0.5)), T.((1.0, 1.0)))

    p = partition(g, BisectFractionPartition(T.((1.0, 0.0)), T(0.2)))
    p1, p2 = p[1], p[2]
    @test nelements(p1) == 20
    @test nelements(p2) == 80

    # all points in p1 are to the left of p2
    pts1 = [centroid(p1, i) for i in 1:nelements(p1)]
    pts2 = [centroid(p2, i) for i in 1:nelements(p2)]
    X1 = reduce(hcat, coordinates.(pts1))
    X2 = reduce(hcat, coordinates.(pts2))
    M1 = maximum(X1, dims=2)
    m2 = minimum(X2, dims=2)
    @test all(X1[1, j] < m2[1] for j in 1:size(X1, 2))
    @test all(X2[1, j] > M1[1] for j in 1:size(X2, 2))

    # flipping normal direction is equivalent to swapping subsets
    p₁ = partition(g, BisectFractionPartition(T.((1.0, 0.0)), T(0.2)))
    p₂ = partition(g, BisectFractionPartition(T.((-1.0, 0.0)), T(0.8)))
    @test nelements(p₁[1]) == nelements(p₂[2]) == 20
    @test nelements(p₁[2]) == nelements(p₂[1]) == 80

    # reproducible results with rng
    rng = StableRNG(123)
    g = CartesianGrid{T}(10, 10)
    p1 = partition(rng, g, BisectFractionPartition(T.((1, 0)), T(0.5)))
    rng = StableRNG(123)
    p2 = partition(rng, g, BisectFractionPartition(T.((1, 0)), T(0.5)))
    @test p1 == p2
  end

  @testset "BallPartition" begin
    pset = PointSet(T[
      0 1 1 0 0.2
      0 0 1 1 0.2
    ])

    # 3 balls with 1 point, and 1 ball with 2 points
    p = partition(pset, BallPartition(T(0.5)))
    n = nelements.(p)
    @test length(p) == 4
    @test count(i -> i == 1, n) == 3
    @test count(i -> i == 2, n) == 1
    @test setify(indices(p)) == setify([[1, 5], [2], [3], [4]])

    # 5 balls with 1 point each
    p = partition(pset, BallPartition(T(0.2)))
    @test length(p) == 5
    @test all(nelements.(p) .== 1)
    @test setify(indices(p)) == setify([[1], [2], [3], [4], [5]])

    # reproducible results with rng
    rng = StableRNG(123)
    g = CartesianGrid{T}(10, 10)
    p1 = partition(rng, g, BallPartition(T(2)))
    rng = StableRNG(123)
    p2 = partition(rng, g, BallPartition(T(2)))
    @test p1 == p2
  end

  @testset "PlanePartition" begin
    g = CartesianGrid((3, 3), T.((-0.5, -0.5)), T.((1.0, 1.0)))
    p = partition(g, PlanePartition(T.((0, 1))))
    @test setify(indices(p)) == setify([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

    g = CartesianGrid((4, 4), T.((-0.5, -0.5)), T.((1.0, 1.0)))
    p = partition(g, PlanePartition(T.((0, 1))))
    @test setify(indices(p)) == setify([1:4, 5:8, 9:12, 13:16])

    # reproducible results with rng
    rng = StableRNG(123)
    g = CartesianGrid{T}(10, 10)
    p1 = partition(rng, g, PlanePartition(T.((1, 0))))
    rng = StableRNG(123)
    p2 = partition(rng, g, PlanePartition(T.((1, 0))))
    @test p1 == p2
  end

  @testset "PredicatePartition" begin
    g = CartesianGrid((3, 3), T.((-0.5, -0.5)), T.((1.0, 1.0)))

    # partition even from odd locations
    pred(i, j) = iseven(i + j)
    partitioner = PredicatePartition(pred)
    p = partition(g, partitioner)
    @test setify(indices(p)) == setify([1:2:9, 2:2:8])

    # reproducible results with rng
    rng = StableRNG(123)
    g = CartesianGrid{T}(10, 10)
    p1 = partition(rng, g, partitioner)
    rng = StableRNG(123)
    p2 = partition(rng, g, partitioner)
    @test p1 == p2
  end

  @testset "SpatialPredicatePartition" begin
    g = CartesianGrid((10, 10), T.((-0.5, -0.5)), T.((1.0, 1.0)))

    # check if there are 100 partitions, each one having only 1 point
    sp = SpatialPredicatePartition((x, y) -> norm(x - y) < T(1))
    s = indices(partition(g, sp))
    @test length(s) == 100
    nelms = [nelements(d) for d in partition(g, sp)]
    @test all(nelms .== 1)

    # defining a predicate to check if points x and y belong to the square [0.,5.]x[0.,5.]
    pred(x, y) = all(T[0, 0] .<= x .<= T[5, 5]) && all(T[0, 0] .<= y .<= T[5, 5])
    sp = SpatialPredicatePartition(pred)
    p = partition(g, sp)
    s = indices(p)
    n = nelements.(p)

    # There will be 65 subsets:
    # 1 subset with 36 points (inside square [0.,5.]x[0.,5.])
    # 64 subsets with only 1 point inside each of them
    @test length(s) == 65
    @test maximum(length.(s)) == 36
    @test count(i -> i == 1, n) == 64
    @test count(i -> i == 36, n) == 1

    # reproducible results with rng
    rng = StableRNG(123)
    g = CartesianGrid{T}(10, 10)
    p1 = partition(rng, g, sp)
    rng = StableRNG(123)
    p2 = partition(rng, g, sp)
    @test p1 == p2
  end

  @testset "ProductPartition" begin
    g = CartesianGrid((100, 100), T.((-0.5, -0.5)), T.((1.0, 1.0)))
    bm = BlockPartition(T(10), T(10))
    bn = BlockPartition(T(5), T(5))
    bmn = ProductPartition(bm, bn)

    # Bm x Bn = Bn with m > n
    s1 = indices(partition(g, bmn))
    s2 = indices(partition(g, bn))
    @test setify(s1) == setify(s2)

    # pXp=p (for deterministic p)
    for p in [BlockPartition(T(10), T(10)), BisectFractionPartition(T.((0.1, 0.1)))]
      pp = ProductPartition(p, p)
      s1 = indices(partition(g, pp))
      s2 = indices(partition(g, p))
      @test setify(s1) == setify(s2)
    end

    # reproducible results with rng
    rng = StableRNG(123)
    g = CartesianGrid{T}(10, 10)
    p1 = partition(rng, g, bmn)
    rng = StableRNG(123)
    p2 = partition(rng, g, bmn)
    @test p1 == p2
  end

  @testset "HierarchicalPartition" begin
    g = CartesianGrid((100, 100), T.((-0.5, -0.5)), T.((1.0, 1.0)))
    bm = BlockPartition(T(10), T(10))
    bn = BlockPartition(T(5), T(5))
    bmn = HierarchicalPartition(bm, bn)

    # Bn -> Bm = Bm with m > n
    s1 = indices(partition(g, bmn))
    s2 = indices(partition(g, bn))
    @test setify(s1) == setify(s2)

    # reproducible results with rng
    rng = StableRNG(123)
    g = CartesianGrid{T}(10, 10)
    p1 = partition(rng, g, bmn)
    rng = StableRNG(123)
    p2 = partition(rng, g, bmn)
    @test p1 == p2
  end

  @testset "Mixed Tests" begin
    g = CartesianGrid((100, 100), T.((-0.5, -0.5)), T.((1.0, 1.0)))
    bm = BlockPartition(T(10), T(10))
    bn = BlockPartition(T(5), T(5))
    bmn = ProductPartition(bm, bn)
    hmn = HierarchicalPartition(bm, bn)

    # Bm*Bn = Bm->Bn
    s1 = indices(partition(g, bmn))
    s2 = indices(partition(g, hmn))
    @test setify(s1) == setify(s2)
  end
end
